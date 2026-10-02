"""Lossless compact models must preserve full fusion and error behavior."""

import json
import shutil
import sqlite3
from pathlib import Path

import pytest

from fusion_function import ReferenceDatabase, annotate_fusion_domains, data
from .preprocessing_fixture import fixture

FIXTURE = Path(__file__).parent / "data/reference.sqlite"
GOLDEN = json.loads((Path(__file__).parent / "data/fusion_golden.json").read_text())


@pytest.fixture
def compressed_reference(tmp_path: Path) -> Path:
    """Encode independent regression models without changing their contents."""
    path = tmp_path / "reference.sqlite"
    shutil.copyfile(FIXTURE, path)
    with sqlite3.connect(path) as db:
        rows = list(db.execute("SELECT * FROM ff_transcripts"))
        db.execute("DROP TABLE ff_transcripts")
        db.execute(
            "CREATE TABLE ff_transcripts (transcript_id TEXT NOT NULL PRIMARY KEY, "
            "version INTEGER, translation_id TEXT, translation_version INTEGER, "
            "status TEXT NOT NULL, payload BLOB NOT NULL)"
        )
        db.executemany(
            "INSERT INTO ff_transcripts VALUES (?,?,?,?,?,?)",
            [(*row[:-1], data.encode_transcript_payload(json.loads(row[-1]))) for row in rows],
        )
        db.execute(
            "INSERT INTO build_metadata VALUES ('transcript_payload_codec', ?)",
            (data.TRANSCRIPT_PAYLOAD_CODEC,),
        )
    return path


@pytest.mark.parametrize("case", GOLDEN, ids=[f"fusion-{i + 1}" for i in range(len(GOLDEN))])
def test_complete_fusion_results_from_compressed_models(
    case: dict, compressed_reference: Path
) -> None:
    """Compare every returned field against independently captured regression outputs."""
    with ReferenceDatabase(compressed_reference) as reference:
        assert (
            annotate_fusion_domains(*case["args"], **case["kwargs"], reference=reference)
            == case["result"]
        )


def test_preprocessing_uses_compact_row_storage_and_readable_errors(tmp_path: Path) -> None:
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries)
    # A rowid lookup confirms this is an ordinary table, regardless of SQL spelling.
    assert db.execute("SELECT rowid FROM ff_transcripts LIMIT 1").fetchone() is not None
    storage = dict(db.execute("SELECT status,typeof(payload) FROM ff_transcripts GROUP BY status"))
    assert storage == {"ready": "blob", "noncoding": "blob", "error": "text"}
    encoded = db.execute(
        "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000001'"
    ).fetchone()[0]
    decoded = data.decode_transcript_payload(encoded)
    assert len(encoded) < len(json.dumps(decoded, separators=(",", ":")).encode()) / 2
    error = json.loads(
        db.execute("SELECT payload FROM ff_transcripts WHERE status='error'").fetchone()[0]
    )
    assert "CDS/peptide length mismatch" in error["error"]
    assert (
        dict(db.execute("SELECT * FROM build_metadata"))["transcript_payload_codec"]
        == "zlib-json-v1"
    )
    assert db.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    db.close()
    with ReferenceDatabase(path) as reference:
        assert reference.get_transcript("ENST00000000001", include_sequence=False) == decoded
        assert reference.get_transcript("ENST00000000004") == error


def test_reject_unknown_transcript_codec(compressed_reference: Path) -> None:
    with sqlite3.connect(compressed_reference) as db:
        db.execute("UPDATE build_metadata SET value='unknown' WHERE key='transcript_payload_codec'")
    with pytest.raises(ValueError, match="Unsupported transcript payload codec"):
        ReferenceDatabase(compressed_reference)
