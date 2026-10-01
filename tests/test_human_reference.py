"""Run automatically when a full reference exists in the configured local cache.

Use --human-reference-db=/path/ensembl.sqlite to select another reference.
These tests use the production-sized FTP reference in read-only mode.
No substitute database, downloads, REST calls or synthetic transcript models.
"""

import json
import re

import pytest

from fusion_function import annotate_fusion_domains

pytestmark = [pytest.mark.integration, pytest.mark.human_reference]

GENES = [
    ("BCR", "ENSG00000186716", "22", 1),
    ("ABL1", "ENSG00000097007", "9", 1),
    ("TP53", "ENSG00000141510", "17", -1),
]


def canonical_transcript(reference, gene_id):
    row = reference.db.execute(
        "SELECT t.stable_id, t.transcript_id FROM ensembl_gene g "
        "JOIN ensembl_transcript t ON t.transcript_id=g.canonical_transcript_id "
        "JOIN ensembl_seq_region r ON r.seq_region_id=t.seq_region_id "
        "JOIN ensembl_coord_system c USING (coord_system_id) "
        "WHERE g.stable_id=? AND g.is_current=1 AND c.version='GRCh38'",
        (gene_id,),
    ).fetchone()
    assert row is not None, f"No canonical transcript for {gene_id}"
    payload = reference.get_transcript(row[0])
    assert "error" not in payload, f"{row[0]}: {payload.get('error')}"
    return row[0], row[1], payload


def test_complete_reference_metadata_and_counts(human_reference_db):
    reference = human_reference_db
    assert reference.metadata["species"] == "homo_sapiens"
    assert reference.metadata["assembly"] == "GRCh38"
    assert re.fullmatch(r"[1-9]\d*", reference.metadata["release"])
    counts = json.loads(reference.metadata["counts"])
    for table in json.loads(reference.metadata["tables"]):
        actual = reference.db.execute(f"SELECT COUNT(*) FROM ensembl_{table}").fetchone()[0]
        assert actual == counts[table], f"Imported count differs for {table}"
    for kind in ("dna", "pep"):
        actual = reference.db.execute(
            "SELECT COUNT(*) FROM sequences WHERE kind=?", (kind,)
        ).fetchone()[0]
        assert actual == counts[kind + "_sequences"] > 0
    actual = dict(
        reference.db.execute("SELECT status, COUNT(*) FROM ff_transcripts GROUP BY status")
    )
    prepared = json.loads(reference.metadata["preprocessing_counts"])
    for status in ("ready", "error", "noncoding"):
        assert actual.get(status, 0) == prepared[status]
    primary = {str(number) for number in range(1, 23)} | {"X", "Y", "MT"}
    stored = {
        row[0] for row in reference.db.execute("SELECT sequence_id FROM sequences WHERE kind='dna'")
    }
    assert primary <= stored


@pytest.mark.parametrize(
    "chromosome,length",
    [("1", 248956422), ("9", 138394717), ("10", 133797422), ("22", 50818468), ("MT", 16569)],
)
def test_known_grch38_chromosome_lengths(human_reference_db, chromosome, length):
    row = human_reference_db.db.execute(
        "SELECT length FROM sequences WHERE kind='dna' AND sequence_id=?", (chromosome,)
    ).fetchone()
    assert row == (length,)


def test_all_dna_records_match_grch38_core(human_reference_db):
    unmatched = human_reference_db.db.execute(
        "SELECT s.sequence_id, s.length FROM sequences s WHERE s.kind='dna' AND NOT EXISTS ("
        "SELECT 1 FROM ensembl_seq_region r JOIN ensembl_coord_system c USING (coord_system_id) "
        "WHERE c.version='GRCh38' AND r.name=s.sequence_id AND r.length=s.length) LIMIT 10"
    ).fetchall()
    assert not unmatched, f"DNA records without a matching GRCh38 core region: {unmatched}"


def test_reported_hschr10_cross_assembly_collision(human_reference_db):
    reference = human_reference_db
    row = reference.db.execute(
        "SELECT length FROM sequences WHERE kind='dna' AND sequence_id='HSCHR10_1_CTG2'"
    ).fetchone()
    assert row == (309802,)
    regions = reference.db.execute(
        "SELECT c.version,r.length FROM ensembl_seq_region r "
        "JOIN ensembl_coord_system c USING (coord_system_id) WHERE r.name='HSCHR10_1_CTG2'"
    ).fetchall()
    assert ("GRCh38", 309802) in regions
    if reference.metadata["release"] == "116":
        assert ("GRCh37", 135582047) in regions
    sequence = reference.sequence("dna", "HSCHR10_1_CTG2", 1, 1000)
    assert len(sequence) == 1000
    assert set(sequence) <= set("ACGTRYSWKMBDHVN")


def test_sequence_chunks_and_reverse_complement(human_reference_db):
    # The production chunk size is 1 MiB. Check an interval crossing a boundary.
    boundary = int(human_reference_db.metadata["sequence_chunk_size"])
    start, end = boundary - 50, boundary + 50
    joined = human_reference_db.sequence("dna", "22", start, end)
    separate = human_reference_db.sequence(
        "dna", "22", start, boundary
    ) + human_reference_db.sequence("dna", "22", boundary + 1, end)
    assert joined == separate
    assert len(joined) == end - start + 1
    expected = joined.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
    assert human_reference_db.sequence("dna", "22", start, end, -1) == expected


@pytest.mark.parametrize("gene,gene_id,chromosome,strand", GENES, ids=[row[0] for row in GENES])
def test_real_transcript_exons_cds_and_peptide(
    human_reference_db, gene, gene_id, chromosome, strand
):
    reference = human_reference_db
    transcript_id, internal_id, payload = canonical_transcript(reference, gene_id)
    assert payload["assembly_name"] == "GRCh38"
    assert payload["chromosome"] == chromosome
    assert payload["strand"] == strand
    sequence = payload["premrna_sequence"]
    assert len(sequence) == payload["premrna_length"]
    assert (
        len(sequence) == payload["transcript_genomic_end"] - payload["transcript_genomic_start"] + 1
    )
    forward = reference.sequence(
        "dna", chromosome, payload["transcript_genomic_start"], payload["transcript_genomic_end"]
    )
    expected = (
        forward
        if strand == 1
        else forward.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
    )
    assert sequence == expected
    exons = reference.db.execute(
        "SELECT et.rank,e.seq_region_start,e.seq_region_end,e.seq_region_strand "
        "FROM ensembl_exon_transcript et JOIN ensembl_exon e USING(exon_id) "
        "WHERE et.transcript_id=? ORDER BY et.rank",
        (internal_id,),
    ).fetchall()
    assert [
        (exon["exon_number"], exon["genomic_start"], exon["genomic_end"], strand)
        for exon in payload["transcript_exons"]
    ] == exons
    assert len(payload["splice_sites"]) == 2 * (len(exons) - 1)
    cds_length = sum(block["cds_end"] - block["cds_start"] + 1 for block in payload["cds_blocks"])
    protein_version = reference.db.execute(
        "SELECT tr.version FROM ensembl_translation tr JOIN ensembl_transcript t "
        "ON t.canonical_translation_id=tr.translation_id WHERE t.transcript_id=?",
        (internal_id,),
    ).fetchone()[0]
    peptide = reference.sequence("pep", f"{payload['translation_id']}.{protein_version}")
    assert len(peptide) == payload["protein_length"]
    assert cds_length in (3 * len(peptide), 3 * len(peptide) + 3)
    # Versioned IDs must round-trip on the actual release, not an API fixture.
    version = reference.db.execute(
        "SELECT version FROM ff_transcripts WHERE transcript_id=?", (transcript_id,)
    ).fetchone()[0]
    assert reference.get_transcript(f"{transcript_id}.{version}") == payload
    assert "error" in reference.get_transcript(f"{transcript_id}.{version + 1}")


@pytest.mark.parametrize("gene,gene_id,chromosome,strand", GENES, ids=[row[0] for row in GENES])
def test_real_domains_have_interpro_metadata(human_reference_db, gene, gene_id, chromosome, strand):
    reference = human_reference_db
    _, _, payload = canonical_transcript(reference, gene_id)
    integrated = [feature for feature in payload["protein_features"] if feature["interpro_id"]]
    assert integrated, f"No integrated protein features for {gene}"
    for feature in integrated:
        entry = reference.get_interpro_annotation(feature["interpro_id"])
        assert entry is not None
        assert entry["entry_type"]
        assert feature["interpro_entry_type"] == entry["entry_type"]
        assert feature["interpro_name"] == entry["name"]
        assert 1 <= feature["start"] <= feature["end"] <= payload["protein_length"]


def test_bcr_abl1_annotation_uses_complete_reference(human_reference_db):
    result = annotate_fusion_domains(
        transcript1_id="ENST00000305877",
        transcript2_id="ENST00000318560",
        breakpoint1="22:23289621",
        breakpoint2="9:130854064",
        gene1_terminus="N",
        gene2_terminus="C",
        reference=human_reference_db,
    )
    assert "error" not in result, result
    assert result["frame_status"] in {"in_frame", "out_of_frame", "in_frame/out_of_frame"}
    assert result["domains"]
    assert {domain["transcript_id"] for domain in result["domains"]} <= {
        "ENST00000305877",
        "ENST00000318560",
    }
    assert any(
        domain["transcript_id"] == "ENST00000318560"
        and domain["breakpoint_based_status"] != "excluded"
        for domain in result["domains"]
    )
