"""Current InterPro lists can omit accessions retained by Ensembl."""

import json

import pytest

from fusion_function import data
from .preprocessing_fixture import fixture, uniprot_file


def entry_list(path, rows):
    path.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n" + rows)
    return path


def historical_reference(tmp_path, monkeypatch):
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    db.execute("UPDATE ensembl_interpro SET interpro_ac='IPR001423'")
    db.commit()
    current = entry_list(tmp_path / "current.list", "IPR1\tDomain\tCurrent entry\n")
    newer = entry_list(tmp_path / "110.list", "IPR1\tDomain\tCurrent entry\n")
    archive = entry_list(tmp_path / "109.list", "IPR001423\tConserved_site\tHistorical site\n")
    fetched = []
    monkeypatch.setattr(data, "listing", lambda _: {"109.0", "110.0", "README"})

    def download(url, directory):
        fetched.append((url, directory))
        local = newer if "/110.0/" in url else archive
        return local, data.sha256_file(local)

    monkeypatch.setattr(data, "download", download)
    return path, db, current, archive, fetched


def test_archived_types_are_recovered_and_remain_available_offline(tmp_path, monkeypatch):
    path, db, current, archive, fetched = historical_reference(tmp_path, monkeypatch)
    cache = tmp_path / "metadata"
    counts = data.preprocess_reference(db, current, interpro_archive_cache=cache)
    assert counts["ready"] == 2
    assert [url for url, _ in fetched] == [
        data.INTERPRO_RELEASES_URL + version + "/entry.list" for version in ("110.0", "109.0")
    ]
    assert fetched[-1][1] == cache / "releases" / "109.0"
    metadata = dict(db.execute("SELECT * FROM build_metadata"))
    assert json.loads(metadata["interpro_historical_entries"]) == {"IPR001423": "109.0"}
    source = json.loads(metadata["interpro_historical_sources"])[0]
    assert source["sha256"] == data.sha256_file(archive)
    assert source["url"] == data.INTERPRO_RELEASES_URL + "109.0/entry.list"
    db.close()
    monkeypatch.setattr(
        data, "download", lambda *_: pytest.fail("Offline reprocessing downloaded metadata")
    )
    assert (
        data.main(
            [
                "--preprocess-only",
                str(path),
                "--interpro-entries",
                str(current),
                "--uniprot-features",
                str(uniprot_file(tmp_path)),
            ]
        )
        == 0
    )
    with data.ReferenceReader(path) as reference:
        assert reference.get_interpro_annotation("IPR001423")["entry_type"] == "conserved_site"
        assert (
            reference.metadata["interpro_historical_entries"]
            == metadata["interpro_historical_entries"]
        )
        assert (
            reference.get_transcript("ENST00000000001")["protein_features"][0][
                "interpro_entry_type"
            ]
            == "conserved_site"
        )


def test_current_metadata_takes_precedence_without_archive_lookup(tmp_path, monkeypatch):
    db = fixture(tmp_path / "reference.sqlite")
    current = entry_list(tmp_path / "current.list", "IPR1\tDomain\tCurrent entry\n")
    monkeypatch.setattr(data, "listing", lambda *_: pytest.fail("Unnecessary archive lookup"))
    data.preprocess_reference(db, current, interpro_archive_cache=tmp_path / "metadata")
    assert db.execute("SELECT entry_type FROM ff_interpro WHERE interpro_id='IPR1'").fetchone() == (
        "domain",
    )
    db.close()


def test_unresolved_archive_metadata_fails_before_transcript_loop_and_rolls_back(
    tmp_path, monkeypatch
):
    path, db, current, archive, fetched = historical_reference(tmp_path, monkeypatch)
    data.preprocess_reference(db, archive)
    before = list(db.execute("SELECT * FROM ff_transcripts"))
    before_metadata = list(db.execute("SELECT * FROM build_metadata"))
    incomplete = entry_list(tmp_path / "incomplete.list", "IPR001423\t\tIncomplete\n")
    entry_list(archive, "IPR1\tDomain\tUnrelated\n")
    monkeypatch.setattr(
        data, "reference_structure", lambda *_: pytest.fail("Transcript loop started")
    )
    with pytest.raises(ValueError, match="transcript preprocessing has not started"):
        data.preprocess_reference(db, incomplete, interpro_archive_cache=tmp_path / "metadata")
    assert list(db.execute("SELECT * FROM ff_transcripts")) == before
    assert list(db.execute("SELECT * FROM build_metadata")) == before_metadata
    db.close()
