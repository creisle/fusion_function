import gzip
import io
import logging
import sqlite3
import threading

import pytest

from fusion_function import data
from .preprocessing_fixture import fixture, uniprot_file


@pytest.fixture
def bars(monkeypatch):
    """Exercise real tqdm counters and cleanup, even under pytest capture."""
    original = data.progress
    created = []

    def recording(*args, **kwargs):
        kwargs.update(disable=False, file=io.StringIO())
        bar = original(*args, **kwargs)
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


def test_sqlite_activity_has_no_invented_total_and_resets_handler(bars):
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
    assert bars[-1].total is None
    assert db.handlers[-1] is None
    with pytest.raises(sqlite3.OperationalError):
        data.sqlite_phase(db, "SELECT * FROM missing", "Failing query")
    assert db.handlers[-1] is None
    db.db.close()


def test_sqlite_timer_refreshes_during_blocked_execute_and_stops(monkeypatch, bars):
    """Elapsed output must update even if SQLite never calls a progress hook."""
    monkeypatch.setattr(data, "SQLITE_REFRESH_SECONDS", 0.01)
    refreshed = threading.Event()
    caller = threading.get_ident()
    original = data.progress

    def recording(*args, **kwargs):
        bar = original(*args, **kwargs)
        refresh = bar.refresh

        def tracked_refresh(*args, **kwargs):
            if threading.get_ident() != caller:
                refreshed.set()
            return refresh(*args, **kwargs)

        bar.refresh = tracked_refresh
        return bar

    monkeypatch.setattr(data, "progress", recording)

    class Connection:
        def set_progress_handler(self, callback, interval):
            # Deliberately never invoke it: the independent timer must update.
            pass

        def execute(self, sql):
            assert threading.get_ident() == caller
            assert refreshed.wait(2), "Timer did not refresh while execute was blocked"
            return self

        def fetchall(self):
            return [("ok",)]

    assert data.sqlite_phase(Connection(), "PRAGMA integrity_check", "Blocked check") == [("ok",)]
    assert refreshed.is_set()
    assert bars[-1].total is None
    assert not any(
        thread.name == "fusion-function-sqlite-timer" for thread in threading.enumerate()
    )
    with pytest.raises(RuntimeError, match="failed work"):
        with data.sqlite_activity("Failing work"):
            raise RuntimeError("failed work")
    assert not any(
        thread.name == "fusion-function-sqlite-timer" for thread in threading.enumerate()
    )


def test_redirected_sqlite_activity_logs_heartbeat_and_phase_outcome(monkeypatch, caplog):
    monkeypatch.setattr(data, "SQLITE_REFRESH_SECONDS", 0.01)
    monkeypatch.setattr(data, "SQLITE_LOG_SECONDS", 0.01)
    heartbeat = threading.Event()

    class Capture(logging.Handler):
        def emit(self, record):
            if "still running" in record.getMessage():
                heartbeat.set()

    handler = Capture()
    data.LOG.addHandler(handler)
    token = data.STEP_PREFIX.set("[3/3] ")
    try:
        with caplog.at_level(logging.INFO):
            with data.sqlite_activity("Prepare ordered protein feature stream"):
                assert heartbeat.wait(2), "Redirected output did not receive a heartbeat"
    finally:
        data.STEP_PREFIX.reset(token)
        data.LOG.removeHandler(handler)
    assert "[3/3] Prepare ordered protein feature stream: still running" in caplog.text
    assert "[3/3] Prepare ordered protein feature stream: done" in caplog.text
    assert not any(
        thread.name == "fusion-function-sqlite-timer" for thread in threading.enumerate()
    )


def test_preprocessing_labels_preparation_indexing_and_commit(tmp_path, caplog):
    db = fixture(tmp_path / "reference.sqlite")
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    with caplog.at_level(logging.INFO):
        data.preprocess_reference(db, entries)
    for phase in (
        "Preserve existing InterPro metadata",
        "Replace previous derived transcript and InterPro tables",
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
        with steps.step("Download inputs"):
            with data.progress(total=1, desc="Download") as bar:
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
