"""Controlled sequence fixtures for isolated-gene splicing and initiation."""

from copy import deepcopy

import pytest

from fusion_function import annotate_fusion_domains


TRANSCRIPT = "ENST_TEST"
FEATURES = [("early", 2, 3), ("middle", 4, 6), ("late", 7, 10)]


class PreparedReference:
    """Supply a prepared transcript model; exercise the real annotation pipeline."""

    def __init__(self, transcript):
        self.transcript = transcript

    def get_transcript(self, transcript_id):
        assert transcript_id == TRANSCRIPT
        return deepcopy(self.transcript)


def prepared_reference(cds, strand, *, intron=None, utr5="CCC", features=FEATURES):
    """Build a small transcript with known domains and optional two-exon splicing."""
    assert cds.startswith("ATG") and len(cds) % 3 == 0
    if intron is None:
        sequence = utr5 + cds + "CCC"
        blocks = [(len(utr5) + 1, len(utr5) + len(cds), 1, len(cds))]
    else:
        sequence = utr5 + cds[:12] + intron + cds[12:] + "CCC"
        blocks = [
            (len(utr5) + 1, len(utr5) + 12, 1, 12),
            (len(utr5) + 13 + len(intron), len(utr5) + len(cds) + len(intron), 13, len(cds)),
        ]
    genomic_end = 99 + len(sequence)

    def genomic(position):
        return 99 + position if strand == 1 else genomic_end - position + 1

    cds_blocks = [
        {
            "premrna_start": start,
            "premrna_end": end,
            "cds_start": cds_start,
            "cds_end": cds_end,
            "genomic_start": min(genomic(start), genomic(end)),
            "genomic_end": max(genomic(start), genomic(end)),
        }
        for start, end, cds_start, cds_end in blocks
    ]
    sites = []
    if intron is not None:
        for site_type, position, exon_number in [
            ("donor", blocks[0][1], 1),
            ("acceptor", blocks[1][0], 2),
        ]:
            sites.append(
                {
                    "type": site_type,
                    "premrna_position": position,
                    "genomic_position": genomic(position),
                    "exon_number": exon_number,
                    "disruption_start": genomic(position),
                    "disruption_end": genomic(position),
                }
            )
    protein_features = [
        {
            "feature_id": name,
            "source": "Pfam",
            "interpro_id": f"IPR_TEST_{name}",
            "interpro_entry_type": "domain",
            "interpro_name": name,
            "description": name,
            "start": start,
            "end": end,
            "cds_start": 3 * (start - 1) + 1,
            "cds_end": 3 * end,
        }
        for name, start, end in features
    ]
    transcript = {
        "chromosome": "1",
        "strand": strand,
        "transcript_genomic_start": 100,
        "transcript_genomic_end": genomic_end,
        "premrna_length": len(sequence),
        "premrna_sequence": sequence,
        "cds_blocks": cds_blocks,
        "splice_sites": sites,
        "protein_features": protein_features,
    }
    return PreparedReference(transcript), genomic


def annotate(reference, genomic_position, terminus):
    slot = 1 if terminus == "N" else 2
    result = annotate_fusion_domains(
        **{
            f"transcript{slot}_id": TRANSCRIPT,
            f"breakpoint{slot}": f"1:{genomic_position}",
            f"gene{slot}_terminus": terminus,
            "reference": reference,
        }
    )
    assert "error" not in result, result
    assert result["assumed_disruptive"] is True
    return result, {domain["name"]: domain for domain in result["domains"]}


def start_prediction(result):
    """Get the initiation summary for the single source transcript."""
    assert set(result["translation_start"]) == {TRANSCRIPT}
    return result["translation_start"][TRANSCRIPT]


@pytest.mark.parametrize("strand", [1, -1])
def test_n_terminal_sequence_is_spliced_and_checked_normally(strand):
    cds = "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand, intron="ATGCCC")
    result, domains = annotate(reference, genomic(reference.transcript["premrna_length"]), "N")
    assert result["frame_status"] == "in_frame"
    assert "alternative_start_used" not in result
    assert start_prediction(result) == "native_start_retained"
    assert all(d["post_splicing_status"] == "preserved" for d in domains.values())
    assert all(d["post_translation_status"] == "preserved" for d in domains.values())


@pytest.mark.parametrize("strand", [1, -1])
def test_c_terminal_5utr_keeps_native_start_even_with_upstream_atg(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand, utr5="ATGTAA")
    result, domains = annotate(reference, genomic(1), "C")
    assert result["frame_status"] == "in_frame"
    assert "alternative_start_used" not in result
    assert start_prediction(result) == "native_start_retained"
    assert all(d["post_translation_status"] == "preserved" for d in domains.values())


@pytest.mark.parametrize("strand", [1, -1])
def test_retained_start_coordinate_without_atg_uses_next_start(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA", strand)
    sequence = reference.transcript["premrna_sequence"]
    reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    result, domains = annotate(reference, genomic(1), "C")
    assert start_prediction(result) == "alternative_start_found"
    assert domains["early"]["post_translation_status"] == "translation_start_excluded"
    assert domains["middle"]["post_translation_status"] == "translation_start_disrupted"
    assert domains["late"]["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize(
    "cut_cds", [2, 3, 4], ids=["split-ATG-after-A", "split-ATG-after-AT", "ATG-fully-lost"]
)
def test_next_start_excludes_earlier_domain_and_truncates_overlapping_domain(strand, cut_cds):
    # The next ATG is native residue 5: early residues 2-3 are skipped,
    # middle residues 4-6 lose residue 4, and late residues 7-10 survive.
    reference, genomic = prepared_reference("ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA", strand)
    result, domains = annotate(reference, genomic(3 + cut_cds), "C")
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "in_frame"
    assert all(d["breakpoint_based_status"] == "included" for d in domains.values())
    assert all(d["breakpoint_retained_percent"] == 100.0 for d in domains.values())
    assert all(d["post_splicing_status"] == "preserved" for d in domains.values())
    assert domains["early"]["post_translation_status"] == "translation_start_excluded"
    assert domains["middle"]["post_translation_status"] == "translation_start_disrupted"
    assert domains["late"]["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
def test_alternative_start_search_happens_after_splicing(strand):
    cds = "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand, intron="ATGCCC")
    result, domains = annotate(reference, genomic(7), "C")
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "in_frame"
    assert domains["middle"]["post_splicing_status"] == "preserved"
    assert domains["middle"]["post_translation_status"] == "translation_start_disrupted"
    assert domains["late"]["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
def test_first_atg_is_used_even_when_later_atg_would_preserve_frame(strand):
    # AAT|GAA contains an ATG starting at CDS base 5, before an in-frame ATG
    # at CDS base 13. Do not skip the first candidate to get a favourable result.
    cds = "ATG" + "AAT" + "GAA" + "AAA" + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand)
    result, domains = annotate(reference, genomic(7), "C")
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "out_of_frame"
    assert all(d["post_translation_status"] == "frame_disrupted" for d in domains.values())


@pytest.mark.parametrize("strand", [1, -1])
def test_no_remaining_atg_excludes_translation_of_retained_domains(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    result, domains = annotate(reference, genomic(7), "C")
    assert result["frame_status"] is None
    assert "alternative_start_used" not in result
    assert start_prediction(result) == "alternative_start_not_found"
    assert all(d["breakpoint_based_status"] == "included" for d in domains.values())
    assert all(d["post_splicing_status"] == "preserved" for d in domains.values())
    assert all(
        d["post_translation_status"] == "translation_start_excluded" for d in domains.values()
    )


@pytest.mark.parametrize("strand", [1, -1])
def test_next_start_after_all_domains_excludes_them(strand):
    cds = "ATG" + "AAA" * 8 + "ATG" + "TAA"
    reference, genomic = prepared_reference(
        cds, strand, features=[("early", 2, 3), ("middle", 4, 6), ("late", 7, 9)]
    )
    result, domains = annotate(reference, genomic(7), "C")
    assert result["frame_status"] == "in_frame"
    assert start_prediction(result) == "alternative_start_found"
    assert all(
        d["post_translation_status"] == "translation_start_excluded" for d in domains.values()
    )


@pytest.mark.parametrize("strand", [1, -1])
def test_stop_after_alternative_start_excludes_downstream_domains(strand):
    cds = "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand, intron="ATGTAA")
    # No upstream donor survives this intronic cut, so the retained intronic
    # ATG/TAA remains before the C-terminal coding exon.
    result, domains = annotate(reference, genomic(16), "C")
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "in_frame"
    assert domains["early"]["breakpoint_based_status"] == "excluded"
    assert domains["early"]["post_translation_status"] is None
    assert domains["middle"]["post_translation_status"] == "premature_termination_excluded"
    assert domains["late"]["post_translation_status"] == "premature_termination_excluded"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("terminus", ["N", "C"])
def test_utr_only_retained_portion_reports_start_lost_without_search(strand, terminus):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    position = 2 if terminus == "N" else reference.transcript["premrna_length"] - 1
    result, _ = annotate(reference, genomic(position), terminus)
    assert result["frame_status"] is None
    assert start_prediction(result) == "native_start_lost"


@pytest.mark.parametrize("strand", [1, -1])
def test_n_terminal_lost_start_uses_alternative_start(strand):
    # Splicing must remove the intronic ATG/TAA before selecting the alternative.
    reference, genomic = prepared_reference(
        "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA", strand, intron="ATGTAA"
    )
    sequence = reference.transcript["premrna_sequence"]
    reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    result, domains = annotate(reference, genomic(reference.transcript["premrna_length"]), "N")
    assert result["frame_status"] == "in_frame"
    assert start_prediction(result) == "alternative_start_found"
    assert domains["early"]["post_translation_status"] == "translation_start_excluded"
    assert domains["middle"]["post_translation_status"] == "translation_start_disrupted"
    assert domains["late"]["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
def test_n_terminal_lost_start_reports_unsuccessful_search(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    sequence = reference.transcript["premrna_sequence"]
    reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    result, domains = annotate(reference, genomic(reference.transcript["premrna_length"]), "N")
    assert start_prediction(result) == "alternative_start_not_found"
    assert result["frame_status"] is None
    assert all(
        d["post_translation_status"] == "translation_start_excluded" for d in domains.values()
    )


@pytest.mark.parametrize("strand", [1, -1])
def test_n_terminal_alternative_start_checks_frame(strand):
    # The first alternative ATG starts at source CDS base 5, out of frame.
    cds = "ATG" + "AAT" + "GAA" + "AAA" + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand)
    sequence = reference.transcript["premrna_sequence"]
    reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    result, domains = annotate(reference, genomic(reference.transcript["premrna_length"]), "N")
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "out_of_frame"
    assert all(d["post_translation_status"] == "frame_disrupted" for d in domains.values())


@pytest.mark.parametrize("strand", [1, -1])
def test_n_terminal_cut_within_native_atg_reports_start_lost(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    result, _ = annotate(reference, genomic(5), "N")  # Only AT of the native ATG survives.
    assert result["frame_status"] is None
    assert start_prediction(result) == "alternative_start_not_found"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("termini", [("N", "C"), ("C", "N")])
def test_same_transcript_fusion_has_one_final_product_start(strand, termini):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    positions = {"N": 15, "C": 16}
    result = annotate_fusion_domains(
        TRANSCRIPT,
        TRANSCRIPT,
        f"1:{genomic(positions[termini[0]])}",
        f"1:{genomic(positions[termini[1]])}",
        *termini,
        reference=reference,
    )
    assert result["translation_start"] == {TRANSCRIPT: "native_start_retained"}
    assert result["frame_status"] == "in_frame"


@pytest.mark.parametrize("strand", [1, -1])
def test_same_transcript_fusion_5utr_joins_keep_c_terminal_native_start(strand):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    result = annotate_fusion_domains(
        TRANSCRIPT, TRANSCRIPT, f"1:{genomic(1)}", f"1:{genomic(2)}", "N", "C", reference=reference
    )
    assert start_prediction(result) == "native_start_retained"
    assert result["frame_status"] == "in_frame"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("termini", [("N", "C"), ("C", "N")])
def test_paired_fusion_with_lost_native_start_uses_alternative(strand, termini):
    reference, genomic = prepared_reference("ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA", strand)
    sequence = reference.transcript["premrna_sequence"]
    reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    positions = {"N": 9, "C": 10}
    result = annotate_fusion_domains(
        TRANSCRIPT,
        TRANSCRIPT,
        f"1:{genomic(positions[termini[0]])}",
        f"1:{genomic(positions[termini[1]])}",
        *termini,
        reference=reference,
    )
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "in_frame"
    domains = {domain["name"]: domain for domain in result["domains"]}
    assert domains["early"]["post_translation_status"] == "translation_start_excluded"
    assert domains["middle"]["post_translation_status"] == "translation_start_disrupted"
    assert domains["late"]["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
def test_different_transcript_partners_search_when_both_native_starts_are_lost(strand):
    n_reference, n_genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    c_reference, c_genomic = prepared_reference(
        "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA", strand
    )
    sequence = n_reference.transcript["premrna_sequence"]
    n_reference.transcript["premrna_sequence"] = sequence[:3] + "ACG" + sequence[6:]
    c_transcript = "ENST_OTHER"

    class PairedReference:
        def get_transcript(self, transcript_id):
            return deepcopy(
                {TRANSCRIPT: n_reference.transcript, c_transcript: c_reference.transcript}[
                    transcript_id
                ]
            )

    result = annotate_fusion_domains(
        TRANSCRIPT,
        c_transcript,
        f"1:{n_genomic(9)}",
        f"1:{c_genomic(10)}",
        "N",
        "C",
        reference=PairedReference(),
    )
    assert start_prediction(result) == "alternative_start_found"
    assert result["frame_status"] == "in_frame"
    (late,) = [
        d for d in result["domains"] if d["transcript_id"] == c_transcript and d["name"] == "late"
    ]
    assert late["post_translation_status"] == "preserved"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("with_start", [False, True], ids=["no-ATG", "junction-ATG"])
def test_paired_fusion_searches_across_the_final_junction(strand, with_start):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand)
    sequence = reference.transcript["premrna_sequence"]
    sequence = sequence[:3] + "ACG" + sequence[6:]
    if with_start:
        # The N portion ends with A (position 7), the C portion starts with TG
        # (positions 11-12). Joining them makes an ATG absent from either side.
        sequence = sequence[:10] + "TG" + sequence[12:]
    reference.transcript["premrna_sequence"] = sequence
    result = annotate_fusion_domains(
        TRANSCRIPT, TRANSCRIPT, f"1:{genomic(7)}", f"1:{genomic(11)}", "N", "C", reference=reference
    )
    assert start_prediction(result) == (
        "alternative_start_found" if with_start else "alternative_start_not_found"
    )
    assert result["frame_status"] == ("in_frame" if with_start else "out_of_frame")
    (late,) = [domain for domain in result["domains"] if domain["name"] == "late"]
    assert late["post_translation_status"] == (
        "preserved" if with_start else "translation_start_excluded"
    )


@pytest.mark.parametrize("termini", [("N", "C"), ("C", "N")])
@pytest.mark.parametrize("known_slot", [1, 2])
def test_missing_partner_is_ignored_even_when_it_has_a_terminus(termini, known_slot):
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", 1)
    # Unknown partner metadata cannot turn a single-gene event into N/N or C/C.
    terminus = termini[known_slot - 1]
    result = annotate_fusion_domains(
        **{
            f"transcript{known_slot}_id": TRANSCRIPT,
            f"transcript{3 - known_slot}_id": None,
            f"breakpoint{known_slot}": f"1:{genomic(3 if terminus == 'C' else 35)}",
            "gene1_terminus": terminus,
            "gene2_terminus": terminus,
            "reference": reference,
        }
    )
    assert "error" not in result
    assert result["assumed_disruptive"] is True
    assert start_prediction(result) == "native_start_retained"


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("terminus", ["N", "C"])
def test_alternative_splicing_combines_start_outcomes_with_slashes(strand, terminus):
    # One donor retains the native start; an earlier donor skips it. The
    # downstream coding exon has an ATG, usable in either retained terminus.
    cds = "ATG" + "AAA" * 3 + "ATG" + "AAA" * 5 + "TAA"
    reference, genomic = prepared_reference(cds, strand, intron="CCCCCC")
    for position in (2, 3):
        reference.transcript["splice_sites"].append(
            {
                "type": "donor",
                "premrna_position": position,
                "genomic_position": genomic(position),
                "exon_number": 1,
                "disruption_start": genomic(position),
                "disruption_end": genomic(position),
            }
        )
    position = reference.transcript["premrna_length"] if terminus == "N" else 1
    result, _ = annotate(reference, genomic(position), terminus)
    assert start_prediction(result) == "native_start_retained/alternative_start_found"
    assert "alternative_start_used" not in result


@pytest.mark.parametrize("strand", [1, -1])
def test_alternative_splicing_reports_search_success_and_failure(strand):
    # Alternative acceptors keep or remove the only remaining ATG, located
    # in the intron. Both products lose the native start at the breakpoint.
    reference, genomic = prepared_reference("ATG" + "AAA" * 9 + "TAA", strand, intron="ATGCCC")
    reference.transcript["splice_sites"].append(
        {
            "type": "acceptor",
            "premrna_position": 16,
            "genomic_position": genomic(16),
            "exon_number": 2,
            "disruption_start": genomic(16),
            "disruption_end": genomic(16),
        }
    )
    result, _ = annotate(reference, genomic(7), "C")
    assert start_prediction(result) == "alternative_start_found/alternative_start_not_found"
    # Site input order cannot change the order of the summarized statuses.
    reference.transcript["splice_sites"].reverse()
    reordered, _ = annotate(reference, genomic(7), "C")
    assert reordered["translation_start"] == result["translation_start"]
