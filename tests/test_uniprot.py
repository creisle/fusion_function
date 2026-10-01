"""Sequence/isoform verification and functional-site retention stay offline."""

import gzip
import json
import sqlite3
from pathlib import Path
from typing import cast
from xml.etree import ElementTree as ET

import pytest

from fusion_function import annotate_fusion_domains, data
from fusion_function.reference import ReferenceDatabase
from fusion_function.ensembl import ProteinFeature
from fusion_function.fusion import FusionSide, _annotate_feature_statuses
from fusion_function.uniprot import (
    CuratedProtein,
    NS,
    entry_proteins,
    exact_interval,
    import_features,
)
from .preprocessing_fixture import fixture


AMACR_XML = Path(__file__).parent / "data/AMACR_uniprot.xml.gz"


def small_xml(tmp_path: Path, *, sequence: str = "MKLFIN", fuzzy: bool = False) -> Path:
    """Reviewed human features in the six-residue synthetic reference peptide."""
    path = tmp_path / "sites.xml"
    status = ' status="uncertain"' if fuzzy else ""
    path.write_text(f'''<uniprot xmlns="http://uniprot.org/uniprot">
      <entry dataset="Swiss-Prot"><accession>PTEST</accession>
        <organism><dbReference type="NCBI Taxonomy" id="9606"/></organism>
        <dbReference type="Ensembl" id="ENST00000000001.3">
          <property type="protein sequence ID" value="ENSP00000000001.2"/>
        </dbReference>
        <feature type="active site" description="Proton donor" evidence="1">
          <location><position position="2"{status}/></location>
        </feature>
        <feature type="binding site" evidence="2">
          <location><begin position="4"/><end position="5"/></location>
          <ligand><name>substrate</name></ligand>
        </feature>
        <feature type="motif" description="Targeting signal" evidence="2">
          <location><position position="6"/></location>
        </feature>
        <evidence key="1" type="ECO:0000250"><source>
          <dbReference type="UniProtKB" id="POTHER"/>
        </source></evidence>
        <evidence key="2" type="ECO:0000269"><source>
          <dbReference type="PubMed" id="123456"/>
        </source></evidence>
        <sequence>{sequence}</sequence>
      </entry></uniprot>''')
    return path


def amacr_proteins() -> list[CuratedProtein]:
    with gzip.open(AMACR_XML) as stream:
        entry = ET.parse(stream).getroot().find(NS + "entry")
    return list(entry_proteins(entry))


def test_real_amacr_isoforms_keep_sites_but_not_changed_binding_or_targeting_regions() -> None:
    proteins = amacr_proteins()
    canonical = next(p for p in proteins if "ENSP00000334424" in p["protein_ids"])
    long_isoform = next(p for p in proteins if "ENSP00000371517" in p["protein_ids"])
    shortened = next(p for p in proteins if "ENSP00000371504" in p["protein_ids"])
    assert len(canonical["sequence"]) == 382
    assert len(long_isoform["sequence"]) == 394
    assert long_isoform["sequence"][:377] == canonical["sequence"][:377]
    active = [f for f in long_isoform["features"] if f["feature_type"] == "active_site"]
    assert [f["start"] for f in active] == [122, 152]
    assert all(
        f["evidence"] == [{"code": "ECO:0000250", "source": "UniProtKB", "id": "O06543"}]
        for f in active
    )
    assert not any(
        f["description"] == "Microbody targeting signal" for f in long_isoform["features"]
    )
    assert [f["start"] for f in shortened["features"] if f["feature_type"] == "active_site"] == [
        122
    ]
    assert [
        (f["start"], f["end"]) for f in canonical["features"] if f["feature_type"] == "binding_site"
    ] == [(36, 36), (55, 58), (121, 126)]


def test_real_amacr_reconstructed_isoforms_match_independently_downloaded_sequences() -> None:
    db = sqlite3.connect(":memory:")
    data.create_sequence_tables(db)
    data.import_fasta(db, "pep", AMACR_XML.with_name("AMACR_uniprot_sequences.fa.gz"))
    counts = import_features(db, AMACR_XML)
    assert counts["matched_proteins"] == 2
    assert counts["sequence_mismatches"] == 0
    assert counts["features"] == 11
    db.close()


def test_paper_amacr_fragment_keeps_catalytic_sites_but_loses_n_terminal_binding_sites() -> None:
    protein = next(p for p in amacr_proteins() if "ENSP00000371517" in p["protein_ids"])
    features = [
        cast(
            ProteinFeature,
            {**feature, "cds_start": (feature["start"] - 1) * 3 + 1, "cds_end": feature["end"] * 3},
        )
        for feature in protein["features"]
    ]
    # The paper's AMACR 84–394 fragment retains complete codons from residue 84.
    side = cast(
        FusionSide,
        {
            "side": "transcript2",
            "transcript_id": "ENST00000382085",
            "coding_segments": [{"source_cds_start": 250, "source_cds_end": 1182}],
            "protein_features": features,
        },
    )
    annotated = _annotate_feature_statuses(side, [])
    assert {
        f["start"]: f["breakpoint_based_status"]
        for f in annotated
        if f["feature_type"] == "active_site"
    } == {122: "included", 152: "included"}
    assert {
        f["start"]: f["breakpoint_based_status"]
        for f in annotated
        if f["feature_type"] == "binding_site"
    } == {36: "excluded", 55: "excluded", 121: "included"}


@pytest.mark.parametrize("compressed", [False, True])
def test_exact_sequence_mapping_evidence_and_genomic_segments(
    tmp_path: Path, compressed: bool
) -> None:
    db = fixture(tmp_path / "reference.sqlite")
    path = small_xml(tmp_path)
    if compressed:
        # Download caches need not have a .gz extension.
        compressed_path = tmp_path / "stream"
        compressed_path.write_bytes(gzip.compress(path.read_bytes()))
        path = compressed_path
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries, uniprot_features=path)
    payload = json.loads(
        db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000001'"
        ).fetchone()[0]
    )
    features = [f for f in payload["protein_features"] if f["source"] == "UniProtKB"]
    assert [(f["feature_type"], f["start"], f["end"]) for f in features] == [
        ("active_site", 2, 2),
        ("binding_site", 4, 5),
        ("motif", 6, 6),
    ]
    assert features[0]["genomic_segments"] == [
        {"chromosome": "1", "start": 106, "end": 108, "strand": 1, "assembly_name": "GRCh38"}
    ]
    assert features[1]["genomic_segments"][0]["start"] == 148
    metadata = dict(db.execute("SELECT * FROM build_metadata"))
    assert metadata["feature_annotation_version"] == "2"
    assert json.loads(metadata["uniprot_features_source"])["sha256"] == data.sha256_file(path)
    before = list(db.execute("SELECT * FROM ff_transcripts"))
    data.preprocess_reference(db, entries)  # Persisted lookup works offline.
    assert list(db.execute("SELECT * FROM ff_transcripts")) == before
    db.close()


@pytest.mark.parametrize("mismatch", ["length", "residue", "protein_id", "mouse", "unreviewed"])
def test_unmatched_sequences_species_and_unreviewed_entries_are_not_propagated(
    tmp_path: Path, mismatch: str
) -> None:
    db = fixture(tmp_path / "reference.sqlite")
    path = small_xml(
        tmp_path, sequence={"length": "MKLFINA", "residue": "MKLFIA"}.get(mismatch, "MKLFIN")
    )
    text = path.read_text()
    if mismatch == "protein_id":
        text = text.replace("ENSP00000000001.2", "ENSP00000009999.2")
    elif mismatch == "mouse":
        text = text.replace('id="9606"', 'id="10090"')
    elif mismatch == "unreviewed":
        text = text.replace('dataset="Swiss-Prot"', 'dataset="TrEMBL"')
    path.write_text(text)
    if mismatch in {"mouse", "unreviewed"}:
        with pytest.raises(ValueError, match="no reviewed human"):
            import_features(db, path)
    else:
        counts = import_features(db, path)
        assert counts["features"] == 0
        assert counts["matched_proteins"] == 0
    db.close()


def test_protein_cross_reference_can_map_older_ensembl_ids_after_sequence_verification(
    tmp_path: Path,
) -> None:
    db = fixture(tmp_path / "reference.sqlite")
    db.execute("ALTER TABLE ensembl_xref ADD COLUMN xref_id INTEGER")
    db.execute("INSERT INTO ensembl_xref VALUES ('PTEST','Protein','PTEST',100)")
    db.execute(
        "CREATE TABLE ensembl_object_xref (xref_id INTEGER, ensembl_id INTEGER, ensembl_object_type TEXT)"
    )
    db.execute("INSERT INTO ensembl_object_xref VALUES (100,2,'Translation')")
    path = small_xml(tmp_path)
    assert import_features(db, path)["matched_proteins"] == 2
    assert (
        db.execute(
            "SELECT COUNT(*) FROM ff_uniprot_features WHERE sequence_id='ENSP00000000002.2'"
        ).fetchone()[0]
        == 3
    )
    path.write_text(path.read_text().replace("MKLFIN", "MKLFIA"))
    assert import_features(db, path)["features"] == 0
    db.close()


@pytest.mark.parametrize("replacement,expected,shift", [("", "MKFIN", -1), ("AA", "MKAAFIN", 1)])
def test_isoform_edits_shift_unchanged_features_and_omit_changed_features(
    tmp_path: Path, replacement: str, expected: str, shift: int
) -> None:
    path = small_xml(tmp_path)
    root = ET.parse(path).getroot()
    entry = root.find(NS + "entry")
    comment = ET.SubElement(entry, NS + "comment", {"type": "alternative products"})
    isoform = ET.SubElement(comment, NS + "isoform")
    ET.SubElement(isoform, NS + "id").text = "PTEST-2"
    ET.SubElement(isoform, NS + "sequence", {"type": "described", "ref": "VSP_TEST"})
    variant = ET.SubElement(entry, NS + "feature", {"type": "splice variant", "id": "VSP_TEST"})
    ET.SubElement(variant, NS + "original").text = "L"
    ET.SubElement(variant, NS + "variation").text = replacement
    location = ET.SubElement(variant, NS + "location")
    ET.SubElement(location, NS + "position", {"position": "3"})
    changed = ET.SubElement(
        entry, NS + "feature", {"type": "domain", "description": "Touches altered residue"}
    )
    location = ET.SubElement(changed, NS + "location")
    ET.SubElement(location, NS + "begin", {"position": "2"})
    ET.SubElement(location, NS + "end", {"position": "4"})
    protein = next(p for p in entry_proteins(entry) if "PTEST-2" in p["accessions"])
    assert protein["sequence"] == expected
    assert [(f["start"], f["end"]) for f in protein["features"]] == [
        (2, 2),
        (4 + shift, 5 + shift),
        (6 + shift, 6 + shift),
    ]
    assert not any(f["feature_type"] == "domain" for f in protein["features"])


def test_isoform_specific_features_use_isoform_coordinates(tmp_path: Path) -> None:
    path = small_xml(tmp_path)
    text = path.read_text().replace(
        "<sequence>MKLFIN</sequence>",
        """
      <comment type="alternative products"><isoform><id>PTEST-2</id>
        <sequence type="described" ref="VSP_TEST"/></isoform></comment>
      <feature type="splice variant" id="VSP_TEST"><original>L</original><variation>AA</variation>
        <location><position position="3"/></location></feature>
      <feature type="active site" description="Isoform-specific site">
        <location sequence="PTEST-2"><position position="7"/></location></feature>
      <sequence>MKLFIN</sequence>""",
    )
    proteins = list(entry_proteins(ET.fromstring(text).find(NS + "entry")))
    canonical = next(p for p in proteins if "PTEST" in p["accessions"])
    alternate = next(p for p in proteins if "PTEST-2" in p["accessions"])
    assert not any(f["description"] == "Isoform-specific site" for f in canonical["features"])
    assert (
        next(f for f in alternate["features"] if f["description"] == "Isoform-specific site")[
            "start"
        ]
        == 7
    )


@pytest.mark.parametrize("status", ["uncertain", "unknown", "less than", "greater than"])
def test_fuzzy_coordinates_are_omitted(tmp_path: Path, status: str) -> None:
    element = ET.fromstring(
        f'<feature xmlns="http://uniprot.org/uniprot"><location><position position="2" status="{status}"/></location></feature>'
    )
    assert exact_interval(element) is None
    db = fixture(tmp_path / "reference.sqlite")
    assert import_features(db, small_xml(tmp_path, fuzzy=True))["features"] == 2
    db.close()


def test_invalid_uniprot_replacement_preserves_previous_derived_tables(tmp_path: Path) -> None:
    db = fixture(tmp_path / "reference.sqlite")
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries, uniprot_features=small_xml(tmp_path))
    before = {
        table: list(db.execute(f"SELECT * FROM {table}"))
        for table in ("ff_transcripts", "ff_uniprot_features", "build_metadata")
    }
    bad = tmp_path / "bad.xml"
    bad.write_text("<uniprot")
    with pytest.raises(ET.ParseError):
        data.preprocess_reference(db, entries, uniprot_features=bad)
    for table, records in before.items():
        assert list(db.execute(f"SELECT * FROM {table}")) == records
    db.close()


def test_runtime_retention_uses_sites_and_preserves_old_group_boundaries(tmp_path: Path) -> None:
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries, uniprot_features=small_xml(tmp_path))
    payload = json.loads(
        db.execute(
            "SELECT payload FROM ff_transcripts WHERE transcript_id='ENST00000000001'"
        ).fetchone()[0]
    )
    # Simulate the overlap groups from a previously built database.
    payload["feature_groups"] = [[0, 1], [2], [3], [4]]
    db.execute(
        "UPDATE ff_transcripts SET payload=? WHERE transcript_id='ENST00000000001'",
        (json.dumps(payload),),
    )
    db.commit()
    db.close()
    result = annotate_fusion_domains(
        "ENST00000000001", None, "1:145", "1:900", "N", "", database=path
    )
    assert "error" not in result
    domains = [f for f in result["domains"] if f["interpro_id"] == "IPR1"]
    assert {(f["start"], f["end"]) for f in domains} == {(2, 5), (3, 5)}
    assert {f["breakpoint_retained_percent"] for f in domains} == {33.3, 50.0}
    active = next(f for f in result["domains"] if f["domain_type"] == "active_site")
    assert active["breakpoint_based_status"] == "included"
    assert active["sources"] == ["UniProtKB"]
    assert active["uniprot_accessions"] == ["PTEST"]
    assert active["evidence"][0]["code"] == "ECO:0000250"
    assert (
        next(f for f in result["domains"] if f["domain_type"] == "binding_site")[
            "breakpoint_based_status"
        ]
        == "excluded"
    )


@pytest.mark.integration
@pytest.mark.human_reference
def test_built_amacr_has_verified_functional_sites(human_reference_db: ReferenceDatabase) -> None:
    if human_reference_db.metadata.get("feature_annotation_version") != "2":
        pytest.skip("Reprocess the built database to add reviewed UniProt features")
    transcript = human_reference_db.get_transcript("ENST00000335606")
    assert "error" not in transcript
    sites = [f for f in transcript["protein_features"] if f.get("uniprot_accession") == "Q9UHK6"]
    assert {f["start"] for f in sites if f["feature_type"] == "active_site"} == {122, 152}
    assert {(f["start"], f["end"]) for f in sites if f["feature_type"] == "binding_site"} == {
        (36, 36),
        (55, 58),
        (121, 126),
    }


@pytest.mark.integration
@pytest.mark.human_reference
def test_built_slc45a2_amacr_evaluates_sites_in_final_product(
    human_reference_db: ReferenceDatabase,
) -> None:
    from .test_fusion import _breakpoint_after_exon

    if human_reference_db.metadata.get("feature_annotation_version") != "2":
        pytest.skip("Reprocess the built database to add reviewed UniProt features")
    slc45a2 = human_reference_db.get_transcript("ENST00000296589")
    amacr = human_reference_db.get_transcript("ENST00000335606")
    assert "error" not in slc45a2 and "error" not in amacr
    result = annotate_fusion_domains(
        "ENST00000296589",
        "ENST00000335606",
        _breakpoint_after_exon(slc45a2, 2),
        _breakpoint_after_exon(amacr, 1),
        "N",
        "C",
        reference=human_reference_db,
    )
    assert "error" not in result
    features = [f for f in result["domains"] if f.get("uniprot_accessions") == ["Q9UHK6"]]
    active = [f for f in features if f["domain_type"] == "active_site"]
    assert {f["start"] for f in active} == {122, 152}
    assert all(f["breakpoint_based_status"] == "included" for f in active)
    assert all(f["post_splicing_status"] == "preserved" for f in active)
    assert all("preserved" in (f["post_translation_status"] or "").split("/") for f in active)
    binding = {f["start"]: f for f in features if f["domain_type"] == "binding_site"}
    assert binding[36]["breakpoint_based_status"] == "excluded"
    assert binding[55]["breakpoint_based_status"] == "excluded"
    assert binding[121]["breakpoint_based_status"] == "included"
    assert {f["domain_type"] for f in result["domains"] if f["interpro_id"] == "IPR050509"} == {
        "family"
    }
