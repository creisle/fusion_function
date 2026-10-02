"""Regression checks for Ensembl's virtual leading peptide-codon bases."""

import hashlib
import zlib
from pathlib import Path
from typing import Literal

import pytest

from fusion_function import ReferenceDatabase, annotate_fusion_domains, data
from .preprocessing_fixture import fixture


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("phase", [1, 2])
@pytest.mark.parametrize("has_stop", [False, True])
@pytest.mark.parametrize("terminus", ["N", "C"])
def test_missing_leading_codon_bases_preserve_mapping_and_start_predictions(
    tmp_path: Path, strand: int, phase: int, has_stop: bool, terminus: Literal["N", "C"]
) -> None:
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    tid = 1 if strand == 1 else 2
    transcript_id = f"ENST{tid:011}"
    # The first residue is unknown because its codon lacks `phase` bases.
    # Residues 2..6 are independently encoded as MKLFN; an optional stop follows.
    coding_sequence = "A" * (3 - phase) + "ATGAAACTTTTTAAC" + ("TAA" if has_stop else "")
    db.execute("UPDATE ensembl_exon SET phase=? WHERE exon_id=?", (phase, tid * 10 + 1))
    db.execute(
        "UPDATE ensembl_translation SET seq_start=1, seq_end=? WHERE translation_id=?",
        (len(coding_sequence) - 12, tid),
    )
    # Include a feature in the partial first residue as well as residues 2..5.
    db.execute(
        "UPDATE ensembl_protein_feature SET seq_start=1,seq_end=1 WHERE protein_feature_id=?",
        (tid * 10 + 2,),
    )
    if strand == 1:
        positions = list(range(100, 112)) + list(range(148, 148 + len(coding_sequence) - 12))
    else:
        positions = list(range(159, 147, -1)) + list(
            range(111, 111 - (len(coding_sequence) - 12), -1)
        )
    genome = list(data.fetch_sequence(db, "dna", "1"))
    for position, nucleotide in zip(positions, coding_sequence, strict=True):
        genome[position - 1] = (
            nucleotide if strand == 1 else nucleotide.translate(str.maketrans("ACGT", "TGCA"))
        )
    for kind, sequence_id, sequence in (
        ("dna", "1", "".join(genome)),
        ("pep", f"ENSP{tid:011}.2", "XMKLFN"),
    ):
        encoded = sequence.encode()
        db.execute(
            "UPDATE sequence_chunks SET data=? WHERE kind=? AND sequence_id=?",
            (zlib.compress(encoded), kind, sequence_id),
        )
        db.execute(
            "UPDATE sequences SET sha256=? WHERE kind=? AND sequence_id=?",
            (hashlib.sha256(encoded).hexdigest(), kind, sequence_id),
        )
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    counts = data.preprocess_reference(db, entries)
    assert counts["ready"] == 2 and counts["error"] == 1
    db.close()
    with ReferenceDatabase(path) as reference:
        transcript = reference.get_transcript(transcript_id)
        assert "error" not in transcript
        assert transcript["cds_start_phase"] == phase
        assert transcript["cds_blocks"][0]["cds_start"] == phase + 1
        assert sum(b["cds_end"] - b["cds_start"] + 1 for b in transcript["cds_blocks"]) == len(
            coding_sequence
        )
        peptide = reference.sequence("pep", f"ENSP{tid:011}.2")
        assert peptide == "XMKLFN"
        # Compare actual genomic positions to an independent per-base map.
        for feature in transcript["protein_features"]:
            actual = []
            for segment in feature["genomic_segments"]:
                region = list(range(segment["start"], segment["end"] + 1))
                actual.extend(region if strand == 1 else reversed(region))
            assert (
                actual
                == positions[
                    max(0, (feature["start"] - 1) * 3 - phase) : feature["end"] * 3 - phase
                ]
            )
        # Retain the whole annotated transcript, with either single-gene terminus.
        genomic = 159 if (strand, terminus) in {(1, "N"), (-1, "C")} else 100
        result = annotate_fusion_domains(
            transcript1_id=transcript_id,
            transcript2_id=None,
            breakpoint1=f"1:{genomic}",
            gene1_terminus=terminus,
            reference=reference,
        )
        assert "error" not in result
        assert result["frame_status"] == "in_frame"
        assert result["translation_start"] == {transcript_id: "alternative_start_found"}
        whole = next(feature for feature in result["domains"] if feature["start"] == 2)
        assert whole["breakpoint_based_status"] == "included"
        assert whole["post_splicing_status"] == "preserved"
        assert whole["post_translation_status"] == "preserved"
        partial = next(feature for feature in result["domains"] if feature["start"] == 1)
        assert partial["breakpoint_retained_percent"] == round(100 * (3 - phase) / 3, 1)
        assert partial["post_translation_status"] == "translation_start_excluded"


@pytest.mark.parametrize("phase", [None, -2, 3])
def test_invalid_start_phase_remains_an_error(tmp_path: Path, phase: int | None) -> None:
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    db.execute("UPDATE ensembl_exon SET phase=? WHERE exon_id=11", (phase,))
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries)
    db.close()
    with ReferenceDatabase(path) as reference:
        assert (
            "Invalid translation start exon phase"
            in reference.get_transcript("ENST00000000001")["error"]
        )


@pytest.mark.parametrize("phase", [1, 2])
def test_start_phase_does_not_allow_unrelated_length_mismatches(tmp_path: Path, phase: int) -> None:
    path = tmp_path / "reference.sqlite"
    db = fixture(path)
    db.execute("UPDATE ensembl_exon SET phase=? WHERE exon_id=11", (phase,))
    db.execute("UPDATE ensembl_translation SET seq_start=1 WHERE translation_id=1")
    # 21 real bases plus phase cannot explain a four-residue peptide.
    db.execute("UPDATE sequences SET length=4 WHERE sequence_id='ENSP00000000001.2'")
    entries = tmp_path / "entry.list"
    entries.write_text("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR1\tDomain\tTest domain\n")
    data.preprocess_reference(db, entries)
    db.close()
    with ReferenceDatabase(path) as reference:
        error = reference.get_transcript("ENST00000000001")["error"]
        assert "CDS/peptide length mismatch" in error and f"start phase {phase}" in error
