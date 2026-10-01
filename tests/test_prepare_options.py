"""Preparation must publish usable references and reject misleading options."""

import gzip
import io
import json
import logging
import re
import sqlite3

import pytest

from fusion_function import data
from .preprocessing_fixture import fixture, uniprot_file


def entries_file(tmp_path):
    path = tmp_path / "entry.list"
    path.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    return path


@pytest.mark.parametrize(
    "options",
    [
        ["--raw-only"],
        ["--tables", "transcript"],
        ["--sequence-types", "dna"],
        ["--core-name", "homo_sapiens_core_116_38"],
        ["--no-fetch-interpro"],
        ["--fetch-interpro"],
    ],
)
def test_unsafe_options_rejected_before_network(options, monkeypatch):
    monkeypatch.setattr(data, "open_url", lambda _: pytest.fail("Network used"))
    with pytest.raises(SystemExit) as exc:
        data.main(options)
    assert exc.value.code == 2


@pytest.mark.parametrize(
    "options",
    [
        ["--release", "116"],
        ["--species", "homo_sapiens"],
        ["--force"],
        ["--output", "other.sqlite"],
        ["--base-url", "https://mirror.example/pub/"],
        ["--no-fetch-interpro"],
        ["--fetch-interpro"],
    ],
)
def test_reprocessing_rejects_inapplicable_options(options, monkeypatch):
    monkeypatch.setattr(data, "open_url", lambda _: pytest.fail("Network used"))
    with pytest.raises(SystemExit) as exc:
        data.main(["--preprocess-only", "reference.sqlite", *options])
    assert exc.value.code == 2


def test_local_interpro_reprocessing_preserves_metadata_and_provenance(tmp_path):
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    entries = entries_file(tmp_path)
    data.preprocess_reference(db, entries)
    before = list(db.execute("SELECT * FROM ff_interpro"))
    digest = dict(db.execute("SELECT * FROM build_metadata"))["interpro_entries_sha256"]
    db.close()
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
    with sqlite3.connect(path) as db:
        assert list(db.execute("SELECT * FROM ff_interpro")) == before
        metadata = dict(db.execute("SELECT * FROM build_metadata"))
        assert metadata["interpro_entries_sha256"] == digest
        assert metadata["interpro_entry_types_loaded"] == "true"
        payload = json.loads(
            db.execute(
                "SELECT payload FROM ff_transcripts WHERE status='ready' LIMIT 1"
            ).fetchone()[0]
        )
        assert payload["protein_features"][0]["interpro_entry_type"] == "domain"


def test_missing_interpro_blocks_preparation_and_rolls_back(tmp_path):
    db = fixture(tmp_path / "reference.sqlite")
    data.preprocess_reference(db, entries_file(tmp_path))
    before = list(db.execute("SELECT * FROM ff_transcripts"))
    metadata = list(db.execute("SELECT * FROM build_metadata"))
    bad = tmp_path / "incomplete.list"
    bad.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\t\tIncomplete\n")
    with pytest.raises(ValueError, match="Missing InterPro entry types"):
        data.preprocess_reference(db, bad)
    assert list(db.execute("SELECT * FROM ff_transcripts")) == before
    assert list(db.execute("SELECT * FROM build_metadata")) == metadata
    assert db.execute("SELECT entry_type FROM ff_interpro").fetchone()[0] == "domain"
    db.close()


@pytest.mark.parametrize(
    "key,value",
    [
        ("species", "mus_musculus"),
        ("assembly", "GRCh37"),
        ("format_version", "2"),
        ("release", ""),
        ("preprocessing_version", "2"),
        ("sequence_chunk_size", "1"),
    ],
)
def test_invalid_source_rejected_before_interpro_download(tmp_path, monkeypatch, key, value):
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    db.execute("INSERT OR REPLACE INTO build_metadata VALUES (?,?)", (key, value))
    db.commit()
    db.close()
    monkeypatch.setattr(data, "download", lambda *_: pytest.fail("Downloaded for invalid database"))
    assert data.main(["--preprocess-only", str(path)]) == 1


def test_interrupted_download_logs_path_and_cleans_partial(tmp_path, monkeypatch, caplog):
    class Interrupted(io.BytesIO):
        headers = {}

        def read(self, *_):
            raise KeyboardInterrupt

    monkeypatch.setattr(data, "open_url", lambda _: Interrupted())
    with caplog.at_level(logging.INFO):
        with pytest.raises(KeyboardInterrupt):
            data.download("https://example.test/table.txt.gz", tmp_path)
    assert "Partial download: " + str(tmp_path) in caplog.text
    assert ".part" in caplog.text
    assert not list(tmp_path.iterdir())


def test_full_build_fixed_inputs_paths_and_atomic_replacement(tmp_path, monkeypatch, caplog):
    """Exercise the actual builder using small Ensembl-shaped FTP files."""
    source = fixture(tmp_path / "source.sqlite")
    for name in ("gene", "external_db", "object_xref"):
        source.execute(f"CREATE TABLE ensembl_{name} (id INTEGER)")
    source.execute("CREATE TABLE ensembl_meta (meta_key TEXT, meta_value TEXT)")
    source.execute("INSERT INTO ensembl_meta VALUES ('schema_version','116')")
    core = "homo_sapiens_core_116_38"
    root = data.BASE_URL + "release-116/"
    core_url = root + f"mysql/{core}/"
    schema = []
    responses = {}
    for name in data.TABLES:
        columns = list(source.execute(f"PRAGMA table_info(ensembl_{name})"))
        schema.append(
            f"CREATE TABLE `{name}` (\n"
            + ",\n".join(
                f"  `{c[1]}` {'int' if c[2] == 'INTEGER' else 'varchar(255)'}" for c in columns
            )
            + "\n);"
        )
        rows = list(source.execute(f"SELECT * FROM ensembl_{name}"))
        text = "".join(
            "\t".join(r"\N" if value is None else str(value) for value in row) + "\n"
            for row in rows
        )
        responses[core_url + name + ".txt.gz"] = gzip.compress(text.encode())
    source.close()
    responses[core_url + core + ".sql.gz"] = gzip.compress("\n".join(schema).encode())
    from fusion_function.uniprot import UNIPROT_HUMAN_URL

    responses[UNIPROT_HUMAN_URL] = gzip.compress(uniprot_file(tmp_path).read_bytes())
    responses[data.INTERPRO_ENTRIES_URL] = entries_file(tmp_path).read_bytes()
    listings = {data.BASE_URL: ["release-115", "release-116"], root + "mysql/": [core]}
    for kind, filename, fasta in [
        ("dna", "Homo_sapiens.GRCh38.dna.toplevel.fa.gz", ">1\n" + "ACGT" * 250 + "\n"),
        (
            "pep",
            "Homo_sapiens.GRCh38.pep.all.fa.gz",
            ">ENSP00000000001.2\nMKLFIN\n>ENSP00000000002.2\nMKLFIN\n>ENSP00000000004.2\nMKLFINS\n",
        ),
    ]:
        url = root + f"fasta/homo_sapiens/{kind}/"
        listings[url] = [filename]
        responses[url + filename] = gzip.compress(fasta.encode())
    for url, names in listings.items():
        responses[url] = "".join(f'<a href="{name}">{name}</a>' for name in names).encode()
    requests = []

    class Response(io.BytesIO):
        def __init__(self, body):
            super().__init__(body)
            self.headers = {"Content-Length": str(len(body))}

    def open_url(url):
        requests.append(url)
        return Response(responses[url])

    monkeypatch.setattr(data, "open_url", open_url)
    cache = tmp_path / "cache"
    with caplog.at_level(logging.INFO):
        assert data.main(["--cache-dir", str(cache)]) == 0
    steps = re.findall(r"\[Step (\d+)/(\d+)\]", caplog.text)
    assert steps == [(str(number), "24") for number in range(1, 25)]
    assert "Completed all 24 steps" in caplog.text
    path = cache / "homo_sapiens/GRCh38/release-116/ensembl.sqlite"
    with data.ReferenceReader(path) as reader:
        assert (
            reader.get_transcript("ENST00000000001")["protein_features"][0]["interpro_entry_type"]
            == "domain"
        )
        assert json.loads(reader.metadata["sequence_types"]) == ["dna", "pep"]
        assert set(json.loads(reader.metadata["tables"])) == set(data.TABLES)
    assert str(cache) in caplog.text
    assert str(path) in caplog.text
    assert str(cache / "metadata/interpro") in caplog.text
    assert "Temporary reference database:" in caplog.text
    assert "Partial download:" in caplog.text
    assert not list(cache.rglob("*.part"))
    assert not list(cache.rglob("*.building.*"))
    assert not any("/cdna/" in url or "/cds/" in url or "/ncrna/" in url for url in requests)
    previous = path.read_bytes()
    imported = []
    original_table_import = data.import_table
    original_import = data.import_fasta
    original_validate = data.validate_dna_regions

    def record_import(db, kind, source, **kwargs):
        imported.append(kind)
        return original_import(db, kind, source, **kwargs)

    monkeypatch.setattr(data, "import_fasta", record_import)
    monkeypatch.setattr(
        data,
        "validate_dna_regions",
        lambda *_: (_ for _ in ()).throw(ValueError("DNA validation failed")),
    )
    assert data.main(["--cache-dir", str(cache), "--release", "116", "--force"]) == 1
    assert imported == ["dna"]  # Fail before spending time importing peptide FASTA.
    assert path.read_bytes() == previous
    monkeypatch.setattr(data, "validate_dna_regions", original_validate)
    original_preprocess = data.preprocess_reference
    monkeypatch.setattr(
        data,
        "preprocess_reference",
        lambda *_, **__: (_ for _ in ()).throw(ValueError("injected failure")),
    )
    assert data.main(["--cache-dir", str(cache), "--release", "116", "--force"]) == 1
    assert path.read_bytes() == previous
    assert not list(cache.rglob("*.building.*"))
    checkpoint = path.with_name(path.name + ".building")
    assert checkpoint.is_file()
    # Fix the metadata failure and repeat the command. All tables and completed
    # FASTA files survive, so none should be imported again.
    monkeypatch.setattr(data, "preprocess_reference", original_preprocess)
    monkeypatch.setattr(
        data, "import_table", lambda *_: pytest.fail("Reimported a completed table")
    )
    monkeypatch.setattr(
        data, "import_fasta", lambda *_, **__: pytest.fail("Reimported completed FASTA")
    )
    assert data.main(["--cache-dir", str(cache), "--release", "116", "--force"]) == 0
    assert not checkpoint.exists()
    with data.ReferenceReader(path) as reader:
        assert "error" not in reader.get_transcript("ENST00000000001")
    # A failure during the final check must not repeat transcript preprocessing.
    monkeypatch.setattr(data, "import_table", original_table_import)
    monkeypatch.setattr(data, "import_fasta", original_import)
    original_phase = data.sqlite_phase

    def fail_integrity(db, statement, description):
        if statement == "PRAGMA integrity_check":
            raise ValueError("injected integrity-check failure")
        return original_phase(db, statement, description)

    monkeypatch.setattr(data, "sqlite_phase", fail_integrity)
    before_integrity = path.read_bytes()
    assert data.main(["--cache-dir", str(cache), "--release", "116", "--force"]) == 1
    assert path.read_bytes() == before_integrity
    with sqlite3.connect(checkpoint) as db:
        assert (
            db.execute(
                "SELECT complete FROM build_checkpoints WHERE stage='preprocessing'"
            ).fetchone()[0]
            == 1
        )
    monkeypatch.setattr(data, "sqlite_phase", original_phase)
    monkeypatch.setattr(
        data,
        "preprocess_reference",
        lambda *_, **__: pytest.fail("Repeated completed preprocessing"),
    )
    monkeypatch.setattr(
        data, "import_table", lambda *_: pytest.fail("Reimported a completed table")
    )
    monkeypatch.setattr(
        data, "import_fasta", lambda *_, **__: pytest.fail("Reimported completed FASTA")
    )
    assert data.main(["--cache-dir", str(cache), "--release", "116", "--force"]) == 0
    assert not checkpoint.exists()
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM source_files").fetchone()[0] == len(data.TABLES) + 5
