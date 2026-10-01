"""FTP analysis labels and bulk subfamilies must preserve functional fallbacks."""

import json
import urllib.error

import pytest

from fusion_function import annotate_fusion_domains, data
from fusion_function.reference import ReferenceDatabase
from .preprocessing_fixture import fixture, uniprot_file


AMACR_ROW = (
    "HUMAN|HGNC=451|UniProtKB=Q9UHK6\tQ9UHK6\tAMACR\tPTHR48228:SF5\t"
    "SUCCINYL-COA--D-CITRAMALATE COA-TRANSFERASE\tALPHA-METHYLACYL-COA RACEMASE\n"
)
LEGACY_KCN_ROW = (
    "HUMAN|HGNC=6280|UniProtKB=O95279\t\tPTHR11003:SF241\t"
    "POTASSIUM CHANNEL, SUBFAMILY K\tPOTASSIUM CHANNEL SUBFAMILY K MEMBER 5\t"
    "cation channel activity#GO:0005261\t\t\t\t\n"
)


def panther_reference(tmp_path):
    """Use real Ensembl analysis labels and the real human classification row."""
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    db.execute("ALTER TABLE ensembl_analysis ADD COLUMN db TEXT")
    db.execute("ALTER TABLE ensembl_analysis ADD COLUMN db_version TEXT")
    db.execute(
        "UPDATE ensembl_analysis SET logic_name='hmmpanther', db='PANTHER', db_version='19.0'"
    )
    db.execute("UPDATE ensembl_protein_feature SET hit_name='PTHR48228'")
    db.execute("DELETE FROM ensembl_interpro")
    db.execute("INSERT INTO ensembl_interpro VALUES ('IPR050509', 'PTHR48228')")
    db.execute("ALTER TABLE ensembl_xref ADD COLUMN xref_id INTEGER")
    db.execute("INSERT INTO ensembl_xref VALUES ('Q9UHK6', 'AMACR', 'AMACR', 100)")
    db.execute(
        "CREATE TABLE ensembl_object_xref (ensembl_id INTEGER, ensembl_object_type TEXT, xref_id INTEGER)"
    )
    # Only one protein has the accession: sharing a family must not enrich others.
    db.execute("INSERT INTO ensembl_object_xref VALUES (1, 'Translation', 100)")
    db.commit()
    entries = tmp_path / "entry.list"
    entries.write_text(
        "ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR050509\tFamily\tCoenzyme A-transferase family III\n"
    )
    classifications = tmp_path / "PTHR19.0_human"
    classifications.write_text(AMACR_ROW)
    return path, db, entries, classifications


def fusion(path):
    return annotate_fusion_domains(
        transcript1_id="ENST00000000001",
        transcript2_id="ENST00000000002",
        breakpoint1="1:145",
        breakpoint2="1:145",
        gene1_terminus="N",
        gene2_terminus="C",
        database=path,
    )


def test_bulk_subfamily_restores_missing_domain_without_changing_coordinates(tmp_path):
    path, db, entries, classifications = panther_reference(tmp_path)
    data.preprocess_reference(db, entries)
    before = json.loads(
        db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000001'"
        ).fetchone()[0]
    )
    assert before["protein_features"][0]["source"] == "PANTHER"
    assert all(feature["domain_type"] == "family" for feature in fusion(path)["domains"])
    data.preprocess_reference(db, entries, panther_classifications=classifications)
    after = json.loads(
        db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000001'"
        ).fetchone()[0]
    )
    assert after["cds_blocks"] == before["cds_blocks"]
    assert (
        after["protein_features"][0]["genomic_segments"]
        == before["protein_features"][0]["genomic_segments"]
    )
    feature = after["protein_features"][0]
    assert feature["panther_subfamily_id"] == "PTHR48228:SF5"
    assert feature["panther_subfamily_description"] == "ALPHA-METHYLACYL-COA RACEMASE"
    other = json.loads(
        db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000002'"
        ).fetchone()[0]
    )
    assert other["protein_features"][0]["panther_subfamily_id"] is None
    domains = fusion(path)["domains"]
    assert len(domains) == 4
    assert all(feature["domain_type"] == "family" for feature in domains)
    assert domains[0]["name"] == "ALPHA-METHYLACYL-COA RACEMASE"
    assert domains[0]["transcript_id"] == "ENST00000000001"
    db.close()


def test_bulk_release_matches_ensembl_and_records_provenance(tmp_path, monkeypatch):
    path, db, entries, classifications = panther_reference(tmp_path)
    fetched = []

    def download(url, directory):
        fetched.append((url, directory))
        return classifications, data.sha256_file(classifications)

    monkeypatch.setattr(data, "download", download)
    cache = tmp_path / "metadata/panther"
    data.preprocess_reference(db, entries, panther_cache=cache)
    expected = (
        data.PANTHER_CLASSIFICATIONS_URL
        + "19.0/PANTHER_Sequence_Classification_files/PTHR19.0_human"
    )
    assert fetched == [(expected, cache / "19.0")]
    source = json.loads(
        dict(db.execute("SELECT * FROM build_metadata"))["panther_classifications_source"]
    )
    assert source == {
        "version": "19.0",
        "url": expected,
        "sha256": data.sha256_file(classifications),
    }
    db.close()


def test_archived_panther_14_1_filename_and_provenance(tmp_path, monkeypatch):
    path, db, entries, classifications = panther_reference(tmp_path)
    db.execute("UPDATE ensembl_analysis SET db_version='14.1'")
    db.execute("UPDATE ensembl_protein_feature SET hit_name='PTHR11003'")
    db.execute("UPDATE ensembl_interpro SET id='PTHR11003'")
    db.execute("UPDATE ensembl_xref SET dbprimary_acc='O95279' WHERE xref_id=100")
    classifications.write_text(LEGACY_KCN_ROW)
    db.commit()
    fetched = []

    def download(url, directory):
        fetched.append(url)
        if url.endswith("_human"):
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
        assert url.endswith("/14.1/PANTHER_Sequence_Classification_files/PTHR14.1_human_")
        return classifications, data.sha256_file(classifications)

    monkeypatch.setattr(data, "download", download)
    data.preprocess_reference(db, entries, panther_cache=tmp_path / "metadata/panther")
    assert len(fetched) == 2
    source = json.loads(
        dict(db.execute("SELECT * FROM build_metadata"))["panther_classifications_source"]
    )
    assert source["version"] == "14.1"
    assert source["url"] == fetched[-1]
    assert source["sha256"] == data.sha256_file(classifications)
    assert fusion(path)["domains"][0]["name"] == "POTASSIUM CHANNEL SUBFAMILY K MEMBER 5"
    db.close()


@pytest.mark.parametrize(
    "row,expected",
    [
        (LEGACY_KCN_ROW, ("O95279", "PTHR11003:SF241", "POTASSIUM CHANNEL SUBFAMILY K MEMBER 5")),
        (AMACR_ROW, ("Q9UHK6", "PTHR48228:SF5", "ALPHA-METHYLACYL-COA RACEMASE")),
    ],
)
def test_panther_column_layouts_use_subfamily_names_not_family_names_or_go(tmp_path, row, expected):
    # An arbitrary filename ensures local overrides get the same format detection.
    path = tmp_path / "classifications.tsv"
    path.write_text(row)
    assert list(data.panther_entry_rows(path)) == [expected]


@pytest.mark.parametrize("row,id_col,name_col", [(LEGACY_KCN_ROW, 2, 4), (AMACR_ROW, 3, 5)])
def test_panther_layout_handles_leading_unassigned_and_unnamed_records(
    tmp_path, row, id_col, name_col
):
    fields = row.rstrip("\n").split("\t")
    unassigned = fields.copy()
    unassigned[id_col] = ""
    unnamed = fields.copy()
    unnamed[name_col] = ""
    path = tmp_path / "classifications.tsv"
    path.write_text("\t".join(unassigned) + "\n" + "\t".join(unnamed) + "\n" + row + row)
    expected = (fields[0].split("UniProtKB=")[1], fields[id_col], fields[name_col])
    assert list(data.panther_entry_rows(path)) == [expected, expected]


@pytest.mark.parametrize("row,id_col", [(LEGACY_KCN_ROW, 2), (AMACR_ROW, 3)])
@pytest.mark.parametrize("first", [False, True])
def test_panther_column_detection_does_not_hide_invalid_subfamilies(tmp_path, row, id_col, first):
    fields = row.rstrip("\n").split("\t")
    fields[id_col] = "PTHR11003:SFbad"
    path = tmp_path / "classifications.tsv"
    path.write_text(("" if first else row) + "\t".join(fields) + "\n")
    with pytest.raises(ValueError, match="Invalid PANTHER subfamily"):
        list(data.panther_entry_rows(path))


def test_cached_archived_panther_filename_is_used_without_network(tmp_path, monkeypatch):
    directory = tmp_path / "14.1"
    directory.mkdir()
    path = directory / "PTHR14.1_human_"
    path.write_text(AMACR_ROW)
    digest = data.sha256_file(path)
    path.with_name(path.name + ".sha256").write_text(digest)
    monkeypatch.setattr(
        data, "open_url", lambda _: pytest.fail("Cached archived metadata used network")
    )
    found, checksum, url = data.download_panther_classifications("14.1", tmp_path)
    assert found == path
    assert checksum == digest
    assert url.endswith("PTHR14.1_human_")


@pytest.mark.parametrize("status", [403, 429, 500])
def test_panther_other_http_errors_do_not_trigger_filename_fallback(tmp_path, monkeypatch, status):
    fetched = []

    def download(url, directory):
        fetched.append(url)
        raise urllib.error.HTTPError(url, status, "Error", {}, None)

    monkeypatch.setattr(data, "download", download)
    with pytest.raises(urllib.error.HTTPError) as error:
        data.download_panther_classifications("14.1", tmp_path)
    assert error.value.code == status
    assert len(fetched) == 1


def test_missing_panther_release_does_not_use_a_different_version(tmp_path, monkeypatch):
    fetched = []

    def download(url, directory):
        fetched.append(url)
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(data, "download", download)
    with pytest.raises(urllib.error.HTTPError):
        data.download_panther_classifications("14.1", tmp_path)
    assert len(fetched) == 2
    assert all("/14.1/" in url for url in fetched)


def test_saved_bulk_metadata_is_reusable_offline(tmp_path, monkeypatch):
    path, db, entries, classifications = panther_reference(tmp_path)
    data.preprocess_reference(db, entries, panther_classifications=classifications)
    before = list(db.execute("SELECT * FROM ff_transcripts"))
    monkeypatch.setattr(
        data, "download", lambda *_: pytest.fail("Offline preprocessing downloaded data")
    )
    data.preprocess_reference(db, entries)
    assert list(db.execute("SELECT * FROM ff_transcripts")) == before
    db.close()


@pytest.mark.parametrize("mismatch", ["family", "gene-only", "ambiguous"])
def test_unmatched_or_ambiguous_assignments_are_not_propagated(tmp_path, mismatch):
    path, db, entries, classifications = panther_reference(tmp_path)
    if mismatch == "family":
        classifications.write_text(AMACR_ROW.replace("PTHR48228:SF5", "PTHR10000:SF1"))
    elif mismatch == "gene-only":
        db.execute("UPDATE ensembl_object_xref SET ensembl_object_type='Gene'")
    else:
        db.execute("INSERT INTO ensembl_xref VALUES ('QOTHER', 'Other protein', 'OTHER', 101)")
        db.execute("INSERT INTO ensembl_object_xref VALUES (1, 'Translation', 101)")
        classifications.write_text(
            AMACR_ROW + AMACR_ROW.replace("Q9UHK6", "QOTHER").replace(":SF5", ":SF6")
        )
    db.commit()
    data.preprocess_reference(db, entries, panther_classifications=classifications)
    assert all(feature["domain_type"] == "family" for feature in fusion(path)["domains"])
    db.close()


def test_invalid_classification_rolls_back_metadata_and_transcripts(tmp_path):
    path, db, entries, classifications = panther_reference(tmp_path)
    data.preprocess_reference(db, entries, panther_classifications=classifications)
    before = {
        table: list(db.execute(f"SELECT * FROM {table}"))
        for table in ("ff_transcripts", "ff_panther", "build_metadata")
    }
    classifications.write_text(AMACR_ROW + "not a human classification\n")
    with pytest.raises(ValueError, match="Invalid human PANTHER"):
        data.preprocess_reference(db, entries, panther_classifications=classifications)
    for table, rows in before.items():
        assert list(db.execute(f"SELECT * FROM {table}")) == rows
    db.close()


def test_cli_can_reprocess_with_local_bulk_metadata(tmp_path):
    path, db, entries, classifications = panther_reference(tmp_path)
    db.close()
    assert (
        data.main(
            [
                "--preprocess-only",
                str(path),
                "--interpro-entries",
                str(entries),
                "--panther-classifications",
                str(classifications),
                "--uniprot-features",
                str(uniprot_file(tmp_path)),
            ]
        )
        == 0
    )
    with ReferenceDatabase(path) as reference:
        feature = reference.get_transcript("ENST00000000001")["protein_features"][0]
        assert feature["source"] == "PANTHER"
        assert feature["panther_subfamily_id"] == "PTHR48228:SF5"
