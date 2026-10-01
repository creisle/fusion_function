import gzip
import json
import sqlite3

import pytest

from fusion_function import data
from .preprocessing_fixture import fixture


def test_preprocessing_strands_errors_and_metadata(tmp_path):
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    counts = data.preprocess_reference(db, entries)
    assert counts == {"ready": 2, "noncoding": 1, "error": 1, "protein_features": 4}
    db.close()
    with data.ReferenceReader(path) as reader:
        plus = reader.get_transcript("ENST00000000001")
        minus = reader.get_transcript("ENST00000000002")
        assert [
            (b["genomic_start"], b["genomic_end"], b["cds_start"], b["cds_end"])
            for b in plus["cds_blocks"]
        ] == [(103, 111, 1, 9), (148, 156, 10, 18)]
        assert [
            (b["genomic_start"], b["genomic_end"], b["cds_start"], b["cds_end"])
            for b in minus["cds_blocks"]
        ] == [(148, 156, 1, 9), (103, 111, 10, 18)]
        assert (
            minus["premrna_sequence"]
            == plus["premrna_sequence"].translate(str.maketrans("ACGT", "TGCA"))[::-1]
        )
        assert plus["feature_groups"] == [[0, 1]]
        assert plus["protein_features"][0]["interpro_entry_type"] == "domain"
        assert "non-coding" in reader.get_transcript("ENST00000000003")["error"]
        assert "mismatch" in reader.get_transcript("ENST00000000004")["error"]
        assert reader.get_interpro_annotation("IPR1")["name"] == "Test domain"


def test_preprocessing_is_atomic_and_reusable_without_downloads(tmp_path):
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    good = tmp_path / "entry.list"
    good.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, good)
    before = list(db.execute("SELECT * FROM ff_transcripts"))
    entries = tmp_path / "bad.tsv"
    entries.write_text("broken header\n")
    with pytest.raises(ValueError, match="InterPro TSV"):
        data.preprocess_reference(db, entries)
    assert list(db.execute("SELECT * FROM ff_transcripts")) == before
    db.close()
    assert data.main(["--preprocess-only", str(path), "--interpro-entries", str(good)]) == 0


def test_preprocess_only_fetches_interpro_by_default(tmp_path, monkeypatch):
    path = tmp_path / "reference.sqlite"
    fixture(path).close()
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    fetched = []

    def download(url, directory):
        fetched.append(url)
        return entries, data.sha256_file(entries)

    monkeypatch.setattr(data, "download", download)
    assert data.main(["--preprocess-only", str(path), "--cache-dir", str(tmp_path)]) == 0
    assert fetched == [data.INTERPRO_ENTRIES_URL]
    with data.ReferenceReader(path) as reader:
        assert reader.get_interpro_annotation("IPR1")["entry_type"] == "domain"


def test_chunk_boundaries_reverse_complement_and_bounded_cache(tmp_path):
    path = tmp_path / "sequence.fa.gz"
    sequence = b"ACGTN" * 500000
    with gzip.open(path, "wb") as stream:
        stream.write(b">1\n" + sequence + b"\n>scaffold.1\nACGTNRYSWKMBDHVN\n")
    db = sqlite3.connect(":memory:")
    data.create_sequence_tables(db)
    data.import_fasta(db, "dna", path)
    a, b = data.CHUNK_SIZE - 10, data.CHUNK_SIZE + 20
    assert data.fetch_sequence(db, "dna", "1", a, b) == sequence[a - 1 : b].decode()
    assert data.fetch_sequence(db, "dna", "scaffold.1", strand=-1) == "NBDHVKMWSRYNACGT"
    from collections import OrderedDict

    cache = OrderedDict()
    data.fetch_sequence(db, "dna", "1", _chunk_cache=cache, _cache_limit=1)
    assert len(cache) == 1
    db.close()


def test_fasta_resume_reuses_complete_records_and_discards_orphan_chunks(tmp_path, monkeypatch):
    path = tmp_path / "dna.fa.gz"
    with gzip.open(path, "wb") as stream:
        stream.write(b">1\nACGT\n>2\nTTTTCCCC\n")
    monkeypatch.setattr(data, "CHUNK_SIZE", 4)
    compress = data.zlib.compress
    calls = 0

    def interrupt(chunk, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise KeyboardInterrupt
        return compress(chunk, **kwargs)

    monkeypatch.setattr(data.zlib, "compress", interrupt)
    database = tmp_path / "partial.sqlite"
    db = sqlite3.connect(database)
    data.create_sequence_tables(db)
    with pytest.raises(KeyboardInterrupt):
        data.import_fasta(db, "dna", path)
    db.close()
    db = sqlite3.connect(database)
    assert list(db.execute("SELECT sequence_id FROM sequences")) == [("1",)]
    # A committed chunk without a finished-record row must never be trusted.
    db.execute(
        "INSERT OR REPLACE INTO sequence_chunks VALUES ('dna','2',1,?)", (compress(b"XXXX"),)
    )
    db.commit()

    def require_reuse(chunk, **kwargs):
        assert chunk != b"ACGT", "Recompressed an already complete sequence"
        return compress(chunk, **kwargs)

    monkeypatch.setattr(data.zlib, "compress", require_reuse)
    assert data.import_fasta(db, "dna", path, resume=True) == 2
    assert data.fetch_sequence(db, "dna", "1") == "ACGT"
    assert data.fetch_sequence(db, "dna", "2") == "TTTTCCCC"
    assert (
        db.execute("SELECT COUNT(*) FROM sequence_chunks WHERE sequence_id='2'").fetchone()[0] == 2
    )
    db.close()


def test_checkpoint_rejects_different_build_configuration(tmp_path):
    db = sqlite3.connect(tmp_path / "checkpoint.sqlite")
    checkpoint = data.BuildCheckpoint(db, {"release": 100, "assembly": "GRCh38"})
    checkpoint.save("table:gene", source={"url": "source", "sha256": "abc", "bytes": 5}, count=12)
    data.BuildCheckpoint(db, {"release": 100, "assembly": "GRCh38"})
    with pytest.raises(ValueError, match="different source/settings"):
        data.BuildCheckpoint(db, {"release": 116, "assembly": "GRCh38"})
    assert checkpoint.read("table:gene")["count"] == 12
    db.close()


def test_build_lock_prevents_concurrent_writers_and_releases_after_failure(tmp_path):
    path = tmp_path / "reference.prepare.lock"
    with pytest.raises(RuntimeError):
        with data.build_lock(path):
            with pytest.raises(ValueError, match="Another build"):
                with data.build_lock(path):
                    pytest.fail("Concurrent builder acquired the same lock")
            raise RuntimeError("injected failure")
    with data.build_lock(path):
        pass
