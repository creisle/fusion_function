import gzip
import io
import logging
import sqlite3

import pytest

from fusion_function import data

from .preprocessing_fixture import fixture, uniprot_file


@pytest.fixture
def bars(monkeypatch):
    """Exercise real tqdm counters and cleanup, even under pytest capture."""
    original = data.progress
    created = []

    def recording(
        iterable=None,
        *,
        total=None,
        desc="",
        unit="it",
        unit_scale=False,
        unit_divisor=1000,
        disable=None,
        file=None,
    ):
        bar = original(
            iterable,
            total=total,
            desc=desc,
            unit=unit,
            unit_scale=unit_scale,
            unit_divisor=unit_divisor,
            disable=False,
            file=io.StringIO(),
        )
        created.append(bar)
        return bar

    monkeypatch.setattr(data, "progress", recording)
    yield created
    assert all(bar not in data.tqdm._instances for bar in created)


@pytest.mark.parametrize("known_size", [True, False])
def test_download_byte_totals_and_unknown_size(tmp_path, monkeypatch, bars, known_size):
    body = b"x" * (4 * data.CHUNK_SIZE + 17)

    class Response(io.BytesIO):
        headers = {"Content-Length": str(len(body))} if known_size else {}

    monkeypatch.setattr(data, "open_url", lambda _: Response(body))
    path, digest = data.download("https://example.test/transcript.txt.gz", tmp_path)
    bar = bars[-1]
    assert bar.total == (len(body) if known_size else None)
    assert bar.n == len(body)
    assert path.read_bytes() == body
    assert data.sha256_file(path) == digest
    assert not list(tmp_path.glob("*.part"))


def test_multimember_fasta_progress_uses_full_compressed_size(tmp_path, bars, caplog):
    path = tmp_path / "sequences.fa.gz"
    path.write_bytes(gzip.compress(b">1\nACGT\n") + gzip.compress(b">2\nTTAA\n"))
    db = sqlite3.connect(":memory:")
    data.create_sequence_tables(db)
    with caplog.at_level(logging.INFO):
        assert data.import_fasta(db, "dna", path) == 2
    assert bars[-1].total == path.stat().st_size
    assert bars[-1].n == path.stat().st_size
    assert data.fetch_sequence(db, "dna", "2") == "TTAA"
    assert "processing " not in caplog.text
    db.close()


def test_table_import_progress_and_unique_indexes(tmp_path, bars):
    path = tmp_path / "gene.txt.gz"
    with gzip.open(path, "wt") as stream:
        stream.write("".join(f"{number}\tENSG{number}\n" for number in range(10001)))
    db = sqlite3.connect(":memory:")
    assert (
        data.import_table(db, "gene", [("gene_id", "INTEGER"), ("stable_id", "TEXT")], path)
        == 10001
    )
    import_bar, index_bar = bars
    assert import_bar.n == import_bar.total == path.stat().st_size
    assert index_bar.n == index_bar.total == 2
    assert db.execute("SELECT COUNT(*) FROM ensembl_gene").fetchone()[0] == 10001
    db.close()


def test_preprocessing_has_exact_transcript_total(tmp_path, bars):
    db = fixture(tmp_path / "reference.sqlite")
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    counts = data.preprocess_reference(db, entries)
    bar = next(bar for bar in bars if bar.desc == "Preprocess transcripts")
    assert bar.n == bar.total == sum(counts[key] for key in ("ready", "error", "noncoding")) == 4
    db.close()


def test_sqlite_activity_creates_no_bar_and_resets_handler(bars):
    class Connection:
        def __init__(self):
            self.db = sqlite3.connect(":memory:")
            self.handlers = []

        def set_progress_handler(self, callback, interval):
            self.handlers.append(callback)
            self.db.set_progress_handler(callback, interval)

        def execute(self, sql):
            return self.db.execute(sql)

    db = Connection()
    assert data.sqlite_phase(db, "PRAGMA integrity_check", "Check database") == [("ok",)]
    assert not bars
    assert db.handlers[-1] is None
    with pytest.raises(sqlite3.OperationalError):
        data.sqlite_phase(db, "SELECT * FROM missing", "Failing query")
    assert db.handlers[-1] is None
    db.db.close()


@pytest.mark.parametrize(
    ("error", "outcome"),
    [(None, "done"), (RuntimeError("failed work"), "failed"), (KeyboardInterrupt(), "interrupted")],
)
def test_sqlite_activity_logs_only_boundaries(bars, caplog, error, outcome):
    token = data.STEP_PREFIX.set("[3/3] ")
    try:
        with caplog.at_level(logging.INFO):
            if error is None:
                with data.sqlite_activity("Create derived tables"):
                    pass
            else:
                with pytest.raises(type(error)):
                    with data.sqlite_activity("Create derived tables"):
                        raise error
    finally:
        data.STEP_PREFIX.reset(token)
    assert [record.getMessage() for record in caplog.records] == [
        "[3/3] Create derived tables",
        f"[3/3] Create derived tables: {outcome}",
    ]
    assert not bars


@pytest.mark.parametrize(
    "description",
    [
        "Drop previous ff_transcripts",
        "Preprocess transcripts",
        "Check database",
        "Download dna FASTA",
        "Import dna FASTA",
        "Import pep FASTA",
        "Import reviewed human UniProt",
        "Import protein_feature",
        "Index protein_feature",
        "Index prepared transcripts",
        "Analyze database",
        "Read UniProt protein cross-references",
        "Read genome and protein sequence lengths",
        "Prepare ordered exon stream",
        "Prepare ordered protein feature stream",
        "Prepare transcript stream",
        "Commit preprocessed reference",
        "Roll back preprocessing transaction",
    ],
)
def test_long_operations_have_one_advance_notice(description, caplog):
    with caplog.at_level(logging.INFO):
        data.long_process_notice(description)
    assert len(caplog.records) == 1
    assert description in caplog.records[0].getMessage()
    assert "can take" in caplog.records[0].getMessage()


def test_preprocessing_labels_preparation_indexing_and_commit(tmp_path, caplog):
    db = fixture(tmp_path / "reference.sqlite")
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    with caplog.at_level(logging.INFO):
        data.preprocess_reference(db, entries)
    for phase in (
        "Preserve existing InterPro metadata",
        "Drop previous ff_transcripts",
        "Drop previous ff_interpro",
        "Create derived transcript and InterPro tables",
        "Load InterPro names from Ensembl cross-references",
        "Read genome and protein sequence lengths",
        "Prepare ordered exon stream",
        "Prepare ordered protein feature stream",
        "Count current transcripts",
        "Prepare transcript stream",
        "Index prepared transcripts",
        "Commit preprocessed reference",
    ):
        assert f"{phase}: done" in caplog.text
    db.close()


def test_redirected_output_disables_terminal_redraws():
    output = io.StringIO()
    with data.progress(total=10, file=output) as bar:
        assert bar.disable
        bar.update(10)
    assert output.getvalue() == ""


def test_numbered_steps_label_bars_and_reset_on_error(bars, caplog):
    with caplog.at_level(logging.INFO):
        steps = data.BuildSteps(2)
        with steps.step("Download inputs"), data.progress(total=1, desc="Download") as bar:
            assert bar.desc == "[1/2] Download"
            bar.update(1)
        with pytest.raises(ValueError, match="interrupted stage"):
            with steps.step("Import inputs"):
                with data.progress(total=1, desc="Import") as bar:
                    assert bar.desc == "[2/2] Import"
                    raise ValueError("interrupted stage")
    assert data.STEP_PREFIX.get() == ""
    assert "[Step 1/2] Download inputs" in caplog.text
    assert "[Step 2/2] Import inputs" in caplog.text


def test_reprocess_only_has_four_numbered_steps(tmp_path, caplog):
    path = tmp_path / "reference.sqlite"
    fixture(path).close()
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    with caplog.at_level(logging.INFO):
        assert (
            data.main(
                [
                    "--preprocess-only",
                    str(path),
                    "--interpro-entries",
                    str(entries),
                    "--uniprot-features",
                    str(uniprot_file(tmp_path)),
                ]
            )
            == 0
        )
    steps = [record.getMessage() for record in caplog.records if "[Step " in record.getMessage()]
    assert steps == [
        "[Step 1/4] Validate existing reference database",
        "[Step 2/4] Prepare InterPro metadata",
        "[Step 3/4] Prepare reviewed human UniProt features",
        "[Step 4/4] Preprocess transcripts, domains and splice sites",
    ]
    assert "Completed all 4 steps" in caplog.text
