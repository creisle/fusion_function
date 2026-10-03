import json
import shutil
import sqlite3
from pathlib import Path

import pytest

from fusion_function import ReferenceDatabase, annotate_fusion_domains
from fusion_function.cli import main
from fusion_function.reference import get_reference, resolve_database

FIXTURE = Path(__file__).parent / "data/reference.sqlite"
GOLDEN = json.loads((Path(__file__).parent / "data/fusion_golden.json").read_text())


@pytest.mark.parametrize("case", GOLDEN, ids=[f"fusion-{i + 1}" for i in range(len(GOLDEN))])
def test_complete_regression_results(case, reference_db):
    assert (
        annotate_fusion_domains(*case["args"], **case["kwargs"], reference=reference_db)
        == case["result"]
    )


def install_reference(cache, release):
    destination = cache / "homo_sapiens/GRCh38" / f"release-{release}" / "ensembl.sqlite"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(FIXTURE, destination)
    with sqlite3.connect(destination) as db:
        db.execute("UPDATE build_metadata SET value=? WHERE key='release'", (str(release),))
    return destination


def test_default_cache_latest_local_release_and_pinning(tmp_path, monkeypatch):
    monkeypatch.delenv("FUSION_FUNCTION_DB", raising=False)
    monkeypatch.setenv("FUSION_FUNCTION_CACHEDIR", str(tmp_path))
    old = install_reference(tmp_path, 99)
    newest = install_reference(tmp_path, 116)
    assert resolve_database() == newest
    assert resolve_database(release=99) == old
    first = get_reference()
    assert get_reference() is first
    assert get_reference(release=99).path == old
    case = GOLDEN[0]
    assert annotate_fusion_domains(*case["args"], **case["kwargs"]) == case["result"]


def test_explicit_database_and_environment_precedence(tmp_path, monkeypatch):
    a = install_reference(tmp_path, 115)
    b = install_reference(tmp_path, 116)
    monkeypatch.setenv("FUSION_FUNCTION_DB", str(a))
    assert resolve_database() == a
    assert resolve_database(b) == b
    with pytest.raises(ValueError, match="Requested release"):
        ReferenceDatabase(a, release=116)
    monkeypatch.setenv("FUSION_FUNCTION_DB", str(b))
    assert get_reference().path == b


def test_missing_db_versioned_ids_and_options(tmp_path, monkeypatch, reference_db):
    monkeypatch.delenv("FUSION_FUNCTION_DB", raising=False)
    monkeypatch.setenv("FUSION_FUNCTION_CACHEDIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="prepare-data"):
        resolve_database()
    transcript_id, version = reference_db.db.execute(
        "SELECT transcript_id, version FROM ff_transcripts LIMIT 1"
    ).fetchone()
    assert "error" not in reference_db.get_transcript(f"{transcript_id}.{version}")
    assert "error" in reference_db.get_transcript(f"{transcript_id}.{version + 1}")
    assert "error" in reference_db.get_transcript("ENST00000000000")
    with pytest.raises(ValueError, match="cannot be combined"):
        get_reference(reference=reference_db, database=FIXTURE)


def test_reject_unprocessed_or_wrong_species(tmp_path):
    path = tmp_path / "reference.sqlite"
    shutil.copyfile(FIXTURE, path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE build_metadata SET value='mus_musculus' WHERE key='species'")
    with pytest.raises(ValueError, match="human"):
        ReferenceDatabase(path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE build_metadata SET value='homo_sapiens' WHERE key='species'")
        db.execute("DELETE FROM build_metadata WHERE key='preprocessing_version'")
    with pytest.raises(ValueError, match="preprocess-only"):
        ReferenceDatabase(path)


@pytest.mark.parametrize("assembly", [None, "GRCh37"])
def test_assembly_metadata_is_required(tmp_path, assembly):
    path = tmp_path / "reference.sqlite"
    shutil.copyfile(FIXTURE, path)
    with sqlite3.connect(path) as db:
        if assembly is None:
            db.execute("DELETE FROM build_metadata WHERE key='assembly'")
        else:
            db.execute("UPDATE build_metadata SET value=? WHERE key='assembly'", (assembly,))
    with pytest.raises(ValueError, match="GRCh38"):
        ReferenceDatabase(path)


def test_missing_entry_type_fails_explicitly(tmp_path):
    path = tmp_path / "incomplete.sqlite"
    shutil.copyfile(FIXTURE, path)
    case = GOLDEN[0]
    transcript_id = case["kwargs"]["transcript1_id"]
    with sqlite3.connect(path) as db:
        encoded = db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id=?", (transcript_id,)
        ).fetchone()[0]
        payload = json.loads(encoded)
        for feature in payload["protein_features"]:
            if feature["interpro_id"]:
                feature["interpro_entry_type"] = None
        db.execute(
            "UPDATE ff_transcripts SET payload=? WHERE transcript_id=?",
            (json.dumps(payload), transcript_id),
        )
    with ReferenceDatabase(path) as reference:
        with pytest.raises(ValueError, match="Missing InterPro entry type"):
            annotate_fusion_domains(*case["args"], **case["kwargs"], reference=reference)


def test_cli_dispatch_and_automatic_interpro(monkeypatch, tmp_path):
    import fusion_function.data as data

    captured = []
    monkeypatch.setattr(
        data, "build", lambda args: captured.append(args) or tmp_path / "fake.sqlite"
    )
    assert main(["prepare-data", "--from-source", "--release", "116"]) == 0
    assert captured[0].release == 116
    assert captured[0].interpro_entries is None  # No local override: load the default metadata.
    assert captured[0].species == "homo_sapiens"
    for option in ("--raw-only", "--no-fetch-interpro", "--fetch-interpro"):
        with pytest.raises(SystemExit) as exc:
            main(["prepare-data", option])
        assert exc.value.code == 2
    (tmp_path / "entry.list").write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n")
    assert main(["prepare-data", "--interpro-entries", str(tmp_path / "entry.list")]) == 0
    assert captured[1].interpro_entries == tmp_path / "entry.list"
