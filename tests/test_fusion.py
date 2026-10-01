from __future__ import annotations

from typing import cast

import pytest
from fusion_function.reference import ReferenceDatabase

from fusion_function.ensembl import (
    EnsemblError,
    TranscriptProteinFeatureResult,
    get_protein_domains,
)
from fusion_function.fusion import (
    AggregatedDomain,
    FusionDomainResult,
    _finalize_domains,
    _find_stop_position,
    annotate_fusion_domains,
)

GENE2TRANSCRIPT = {
    "BCR": "ENST00000305877",
    "ABL1": "ENST00000318560",
    "EML4": "ENST00000318522",
    "ALK": "ENST00000389048",
    "KIF5B": "ENST00000302418",
    "RET": "ENST00000355710",
    "SLC45A2": "ENST00000296589",
    "AMACR": "ENST00000335606",
    "EWSR1": "ENST00000397938",
    "FLI1": "ENST00000527786",
    "ERG": "ENST00000453032",
}


def _assert_success(result: FusionDomainResult | EnsemblError) -> FusionDomainResult:

    assert "error" not in result

    return cast(FusionDomainResult, result)


def _assert_kinase_preserved(result: FusionDomainResult, transcript_id: str) -> None:

    kinase_domains = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == transcript_id
        and feature["domain_type"] == "domain"
        and "kinase" in feature["name"].lower()
    ]

    assert kinase_domains

    assert any(feature["post_translation_status"] == "preserved" for feature in kinase_domains)


def _get_transcript(
    transcript_id: str, reference: ReferenceDatabase
) -> TranscriptProteinFeatureResult:

    result = get_protein_domains(transcript_id, reference=reference)

    assert "error" not in result, f"{transcript_id} returned an error: {result}"

    return cast(TranscriptProteinFeatureResult, result)


def _breakpoint_after_exon(transcript: TranscriptProteinFeatureResult, exon_number: int) -> str:

    exon = transcript["transcript_exons"][exon_number - 1]

    next_exon = transcript["transcript_exons"][exon_number]

    if transcript["strand"] == 1:
        intron_start = exon["genomic_end"] + 1

        intron_end = next_exon["genomic_start"] - 1

    else:
        intron_start = next_exon["genomic_end"] + 1

        intron_end = exon["genomic_start"] - 1

    assert intron_start <= intron_end

    chromosome = exon["chromosome"] or transcript["chromosome"]

    return f"{chromosome}:{(intron_start + intron_end) // 2}"


def test_find_stop_position_finds_in_frame_stop() -> None:

    assert _find_stop_position("ATGAAATAAGGG", 1) == 7


def test_find_stop_position_finds_first_stop() -> None:

    assert _find_stop_position("ATGAAATAATAG", 1) == 7


def test_find_stop_position_ignores_out_of_frame_stop() -> None:

    assert _find_stop_position("ATGCTAACC", 1) is None


def test_find_stop_position_without_translation_start() -> None:

    assert _find_stop_position("ATGTAA", None) is None


def test_finalize_domains_filters_irrelevant_entry_types() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR011009",
            "name": "Protein kinase-like domain superfamily",
            "domain_type": "homologous_superfamily",
            "sources": ["Pfam"],
            "start": 231,
            "end": 498,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR050198",
            "name": "Non-receptor tyrosine kinases",
            "domain_type": "family",
            "sources": ["Pfam"],
            "start": 67,
            "end": 509,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": None,
            "name": "disorder_prediction",
            "domain_type": None,
            "sources": ["Pfam"],
            "start": 518,
            "end": 996,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert {feature["domain_type"] for feature in result} == {
        "domain",
        "family",
        "homologous_superfamily",
    }


def test_finalize_domains_keeps_functionally_relevant_sites() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR017441",
            "name": "Protein kinase, ATP binding site",
            "domain_type": "binding_site",
            "sources": ["Pfam"],
            "start": 248,
            "end": 271,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR008266",
            "name": "Tyrosine-protein kinase, active site",
            "domain_type": "active_site",
            "sources": ["Pfam"],
            "start": 359,
            "end": 371,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR123456",
            "name": "Conserved catalytic motif",
            "domain_type": "conserved_site",
            "sources": ["Pfam"],
            "start": 350,
            "end": 375,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert {feature["domain_type"] for feature in result} == {
        "domain",
        "binding_site",
        "active_site",
        "conserved_site",
    }


def test_finalize_domains_does_not_merge_sites_with_parent_domain() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR017441",
            "name": "Protein kinase, ATP binding site",
            "domain_type": "binding_site",
            "sources": ["Pfam"],
            "start": 248,
            "end": 271,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR008266",
            "name": "Tyrosine-protein kinase, active site",
            "domain_type": "active_site",
            "sources": ["Pfam"],
            "start": 359,
            "end": 371,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 3

    assert {feature["domain_type"] for feature in result} == {
        "domain",
        "binding_site",
        "active_site",
    }


def test_finalize_domains_keeps_distinct_overlapping_annotations() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Prosite_profiles"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR020635",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Smart"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR001245",
            "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam", "PRINTS"],
            "start": 242,
            "end": 492,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 3

    assert result[0]["domain_type"] == "domain"

    assert result[0]["interpro_id"] in {"IPR000719", "IPR020635", "IPR001245"}


def test_finalize_domains_keeps_unresolved_domain_boundaries() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Prosite_profiles"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": None,
            "name": "Catalytic domain of the Protein Tyrosine Kinase, Abelson kinase",
            "domain_type": "domain",
            "sources": ["CDD"],
            "start": 235,
            "end": 497,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 2

    assert {feature["interpro_id"] for feature in result} == {"IPR000719", None}


def test_finalize_domains_keeps_short_member_boundaries() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR001245",
            "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 492,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": None,
            "name": "Protein kinase fingerprint",
            "domain_type": "domain",
            "sources": ["PRINTS"],
            "start": 353,
            "end": 371,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 2


def test_finalize_domains_keeps_separate_domain_occurrences() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR001849",
            "name": "Pleckstrin homology domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 100,
            "end": 200,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR001849",
            "name": "Pleckstrin homology domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 400,
            "end": 500,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 2


def test_finalize_domains_does_not_transitively_chain_domains() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000001",
            "name": "Domain A",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 100,
            "end": 200,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000002",
            "name": "Domain B",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 130,
            "end": 230,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000003",
            "name": "Domain C",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 190,
            "end": 290,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 3


def test_finalize_domains_does_not_merge_different_feature_types() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000001",
            "name": "Functional domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 100,
            "end": 200,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000002",
            "name": "Conserved site",
            "domain_type": "conserved_site",
            "sources": ["Pfam"],
            "start": 100,
            "end": 200,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 2


def test_finalize_domains_combines_breakpoint_statuses() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "disrupted",
            "breakpoint_retained_percent": 50.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 1

    assert result[0]["breakpoint_based_status"] == "included/disrupted"


def test_finalize_domains_combines_splicing_statuses() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "lost",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 1

    assert result[0]["post_splicing_status"] == "preserved/lost"


def test_finalize_domains_combines_sequence_statuses() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "frame_disrupted",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 1

    assert result[0]["post_translation_status"] == "preserved/frame_disrupted"


def test_finalize_domains_deduplicates_existing_slash_statuses() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved/frame_disrupted",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
    ]

    result = _finalize_domains(domains)

    assert result[0]["post_translation_status"] == "preserved/frame_disrupted"


def test_finalize_domains_keeps_domain_and_family_types() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["AMACR"],
            "interpro_id": "IPR003673",
            "name": "CoA-transferase family III",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 3,
            "end": 350,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        },
        {
            "transcript_id": GENE2TRANSCRIPT["AMACR"],
            "interpro_id": "IPR050509",
            "name": "Coenzyme A-transferase family III",
            "domain_type": "family",
            "sources": ["PANTHER"],
            "start": 1,
            "end": 373,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
            "specific_name": "ALPHA-METHYLACYL-COA RACEMASE",
        },
    ]

    result = _finalize_domains(domains)

    assert len(result) == 2

    assert {feature["domain_type"] for feature in result} == {"domain", "family"}

    assert {feature["name"] for feature in result} == {
        "CoA-transferase family III",
        "ALPHA-METHYLACYL-COA RACEMASE",
    }


def test_finalize_domains_keeps_specific_family_when_no_domain_entry_exists() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["AMACR"],
            "interpro_id": "IPR050509",
            "name": "Coenzyme A-transferase family III",
            "domain_type": "family",
            "sources": ["PANTHER"],
            "start": 1,
            "end": 373,
            "breakpoint_based_status": "disrupted",
            "breakpoint_retained_percent": 83.4,
            "post_splicing_status": None,
            "post_translation_status": None,
            "specific_name": "ALPHA-METHYLACYL-COA RACEMASE",
        }
    ]

    result = _finalize_domains(domains)

    assert len(result) == 1

    assert result[0]["domain_type"] == "family"

    assert result[0]["name"] == "ALPHA-METHYLACYL-COA RACEMASE"

    assert result[0]["breakpoint_based_status"] == "disrupted"

    assert result[0]["breakpoint_retained_percent"] == 83.4


def test_final_domain_schema_is_concise() -> None:

    domains: list[AggregatedDomain] = [
        {
            "transcript_id": GENE2TRANSCRIPT["ABL1"],
            "interpro_id": "IPR000719",
            "name": "Protein kinase domain",
            "domain_type": "domain",
            "sources": ["Pfam"],
            "start": 242,
            "end": 493,
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved",
        }
    ]
    result = _finalize_domains(domains)

    assert set(result[0]) == {
        "transcript_id",
        "interpro_id",
        "name",
        "domain_type",
        "sources",
        "feature_ids",
        "start",
        "end",
        "breakpoint_based_status",
        "breakpoint_retained_percent",
        "post_splicing_status",
        "post_translation_status",
    }

    assert result[0]["sources"] == ["Pfam"]

    assert "domain_id" not in result[0]

    assert "description" not in result[0]

    assert "genomic_segments" not in result[0]

    assert result[0]["breakpoint_retained_percent"] == 100.0


@pytest.mark.integration
def test_bcr_exon13_abl1_exon2_e13a2_p210(reference_db: ReferenceDatabase) -> None:
    """BCR exon 13 :: ABL1 exon 2 (e13a2/p210) is in-frame and preserves ABL1 kinase."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23289621",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    _assert_kinase_preserved(result, GENE2TRANSCRIPT["ABL1"])


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_e14a2_p210(reference_db: ReferenceDatabase) -> None:
    """BCR exon 14 :: ABL1 exon 2 (e14a2/p210) is in-frame and preserves ABL1 kinase."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    _assert_kinase_preserved(result, GENE2TRANSCRIPT["ABL1"])


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_uses_functional_abl1_domain_names(
    reference_db: ReferenceDatabase,
) -> None:
    """BCR::ABL1 reports canonical InterPro names for distinct ABL1 functional domains."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    abl1_domains = {
        feature["interpro_id"]: feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["ABL1"]
        and feature["domain_type"] == "domain"
    }

    expected_names = {
        "IPR001452": "SH3 domain",
        "IPR000980": "SH2 domain",
        "IPR001245": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
    }

    for interpro_id, expected_name in expected_names.items():
        assert interpro_id in abl1_domains

        assert abl1_domains[interpro_id]["name"] == expected_name

        assert abl1_domains[interpro_id]["breakpoint_based_status"] == "included"

        assert abl1_domains[interpro_id]["breakpoint_retained_percent"] == 100.0

        assert abl1_domains[interpro_id]["post_translation_status"] == "preserved"


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_preserves_functional_abl1_sites(
    reference_db: ReferenceDatabase,
) -> None:
    """BCR exon 14 :: ABL1 exon 2 preserves the ABL1 kinase domain and catalytic sites."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    abl1_features = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["ABL1"]
    ]

    assert any(
        feature["domain_type"] == "binding_site"
        and "kinase" in feature["name"].lower()
        and feature["post_translation_status"] == "preserved"
        for feature in abl1_features
    )

    assert any(
        feature["domain_type"] == "active_site"
        and "kinase" in feature["name"].lower()
        and feature["post_translation_status"] == "preserved"
        for feature in abl1_features
    )


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_returns_both_partner_domains(
    reference_db: ReferenceDatabase,
) -> None:
    """BCR exon 14 :: ABL1 exon 2 returns relevant features from both fusion partners."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    transcript_ids = {feature["transcript_id"] for feature in result["domains"]}

    assert GENE2TRANSCRIPT["BCR"] in transcript_ids

    assert GENE2TRANSCRIPT["ABL1"] in transcript_ids


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_excluded_domains_have_no_downstream_status(
    reference_db: ReferenceDatabase,
) -> None:
    """BCR exon 14 :: ABL1 exon 2 does not assign downstream status to excluded domains."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    excluded = [
        feature for feature in result["domains"] if feature["breakpoint_based_status"] == "excluded"
    ]

    assert excluded

    assert all(
        feature["post_splicing_status"] is None and feature["post_translation_status"] is None
        for feature in excluded
    )

    assert all(feature["breakpoint_retained_percent"] == 0.0 for feature in excluded)


@pytest.mark.integration
def test_bcr_exon14_abl1_exon2_public_domain_schema(reference_db: ReferenceDatabase) -> None:
    """BCR exon 14 :: ABL1 exon 2 returns only concise functional feature fields."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130854064",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["domains"]

    assert all(
        set(feature)
        >= {
            "transcript_id",
            "interpro_id",
            "name",
            "domain_type",
            "sources",
            "feature_ids",
            "start",
            "end",
            "breakpoint_based_status",
            "breakpoint_retained_percent",
            "post_splicing_status",
            "post_translation_status",
        }
        for feature in result["domains"]
    )

    assert all(
        feature["domain_type"]
        in {
            "domain",
            "family",
            "homologous_superfamily",
            "binding_site",
            "active_site",
            "conserved_site",
            "motif",
        }
        for feature in result["domains"]
    )


@pytest.mark.integration
def test_eml4_exon13_alk_exon20_variant1(reference_db: ReferenceDatabase) -> None:
    """EML4 exon 13 :: ALK exon 20 (variant 1) is in-frame and preserves ALK kinase."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["EML4"],
            transcript2_id=GENE2TRANSCRIPT["ALK"],
            breakpoint1="2:42295516",
            breakpoint2="2:29223528",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    _assert_kinase_preserved(result, GENE2TRANSCRIPT["ALK"])


@pytest.mark.integration
def test_kif5b_exon15_ret_exon12(reference_db: ReferenceDatabase) -> None:
    """KIF5B exon 15 :: RET exon 12 is in-frame and preserves the RET kinase domain."""

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["KIF5B"],
            transcript2_id=GENE2TRANSCRIPT["RET"],
            breakpoint1="10:32028428",
            breakpoint2="10:43116584",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    _assert_kinase_preserved(result, GENE2TRANSCRIPT["RET"])


@pytest.mark.integration
def test_slc45a2_exon2_amacr_exon2_retains_racemase_domain(reference_db: ReferenceDatabase) -> None:
    """SLC45A2 exon 2 :: AMACR exon 2 retains the AMACR racemase domain."""

    slc45a2 = _get_transcript(GENE2TRANSCRIPT["SLC45A2"], reference_db)

    amacr = _get_transcript(GENE2TRANSCRIPT["AMACR"], reference_db)

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["SLC45A2"],
            transcript2_id=GENE2TRANSCRIPT["AMACR"],
            breakpoint1=_breakpoint_after_exon(slc45a2, 2),
            breakpoint2=_breakpoint_after_exon(amacr, 1),
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert "in_frame" in result["frame_status"]

    racemase_domains = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["AMACR"]
        and "racemase" in feature["name"].lower()
    ]

    assert racemase_domains, (
        f"SLC45A2::AMACR has no racemase feature using {reference_db.path}: {result}"
    )

    retained_racemase_domains = [
        feature for feature in racemase_domains if feature["breakpoint_based_status"] != "excluded"
    ]

    assert retained_racemase_domains

    assert any(
        0.0 < feature["breakpoint_retained_percent"] < 100.0
        for feature in retained_racemase_domains
    )

    # Partial domains still need splicing and translation checks. These statuses
    # describe the surviving portion, while the original domain stays disrupted.
    assert all(feature["post_splicing_status"] is not None for feature in retained_racemase_domains)
    assert any(
        feature["post_translation_status"] is not None for feature in retained_racemase_domains
    )

    partial_domains = [
        feature
        for feature in result["domains"]
        if 0.0 < feature["breakpoint_retained_percent"] < 100.0
    ]
    assert all(feature["post_splicing_status"] is not None for feature in partial_domains)


@pytest.mark.integration
def test_breakpoint_inside_abl1_kinase_marks_domain_disrupted(
    reference_db: ReferenceDatabase,
) -> None:
    """

    Test explicit disruption of a functional domain by placing the fusion

    breakpoint within the ABL1 kinase domain.



    BCR breakpoint 22:23290413 is the canonical BCR exon 14 breakpoint used

    by the e14a2 BCR::ABL1 tests above. The ABL1 breakpoint 9:130874943

    (GRCh38) corresponds to ENST00000318560 c.1161 / codon 387, which lies

    well within the ABL1 protein kinase domain (~aa 242-493).



    The ABL1 breakpoint is deliberately synthetic rather than a known

    BCR::ABL1 breakpoint: it was chosen specifically to cut through the

    kinase domain. Because overlapping InterPro/member annotations can have

    slightly different boundaries and are collapsed in the public result, the

    final status may be "disrupted" or a combined status such as

    "included/disrupted"; either must contain a disrupted component.

    """

    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["BCR"],
            transcript2_id=GENE2TRANSCRIPT["ABL1"],
            breakpoint1="22:23290413",
            breakpoint2="9:130874943",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    disrupted_kinase_domains = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["ABL1"]
        and feature["domain_type"] == "domain"
        and "kinase" in feature["name"].lower()
        and "disrupted" in feature["breakpoint_based_status"].split("/")
    ]

    assert disrupted_kinase_domains


@pytest.mark.integration
def test_ewsr1_exon7_fli1_exon6_type1(reference_db: ReferenceDatabase) -> None:
    """
    EWSR1 exon 7 :: FLI1 exon 6 is the classic type I EWSR1::FLI1 fusion.

    The fixed GRCh38 breakpoints were chosen near the middle of the introns
    producing the canonical exon 7::exon 6 junction:
      - 22:29287870 lies in EWSR1 intron 7/8 (22:29287135-29288605)
      - 11:128793694 lies in FLI1 intron 5/6 (11:128782024-128805365)

    This tests the characteristic disruptive architecture of EWSR1::FLI1:
    the EWSR1 RNA-recognition motif (InterPro IPR000504) is excluded, while
    the FLI1 ETS domain (InterPro IPR000418) remains included and preserved.
    """
    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["EWSR1"],
            transcript2_id=GENE2TRANSCRIPT["FLI1"],
            breakpoint1="22:29287870",
            breakpoint2="11:128793694",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    ewsr1_rna_domains = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["EWSR1"]
        and feature["interpro_id"] == "IPR000504"
    ]
    assert ewsr1_rna_domains
    assert any(
        feature["breakpoint_based_status"] == "excluded"
        and feature["breakpoint_retained_percent"] == 0.0
        for feature in ewsr1_rna_domains
    )

    fli1_ets_domains = [
        feature
        for feature in result["domains"]
        if feature["transcript_id"] == GENE2TRANSCRIPT["FLI1"]
        and feature["interpro_id"] == "IPR000418"
    ]
    assert fli1_ets_domains
    assert any(
        feature["breakpoint_based_status"] == "included"
        and feature["breakpoint_retained_percent"] == 100.0
        and feature["post_translation_status"] == "preserved"
        for feature in fli1_ets_domains
    )


@pytest.mark.integration
def test_ewsr1_exon7_fli1_exon5_type2(reference_db: ReferenceDatabase) -> None:
    """
    EWSR1 exon 7 :: FLI1 exon 5 is the classic type II EWSR1::FLI1 fusion.

    The fixed GRCh38 breakpoints were chosen near the middle of the introns
    producing the canonical exon 7::exon 5 junction:
      - 22:29287870 lies in EWSR1 intron 7/8 (22:29287135-29288605)
      - 11:128777471 lies in FLI1 intron 4/5 (11:128772986-128781957)

    This tests a second canonical EWSR1::FLI1 splice configuration. The
    EWSR1 RNA-recognition motif (InterPro IPR000504) should be excluded,
    while the FLI1 ETS domain (InterPro IPR000418) remains included and preserved.
    """
    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["EWSR1"],
            transcript2_id=GENE2TRANSCRIPT["FLI1"],
            breakpoint1="22:29287870",
            breakpoint2="11:128777471",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    assert any(
        feature["transcript_id"] == GENE2TRANSCRIPT["EWSR1"]
        and feature["interpro_id"] == "IPR000504"
        and feature["breakpoint_based_status"] == "excluded"
        and feature["breakpoint_retained_percent"] == 0.0
        for feature in result["domains"]
    )

    assert any(
        feature["transcript_id"] == GENE2TRANSCRIPT["FLI1"]
        and feature["interpro_id"] == "IPR000418"
        and feature["breakpoint_based_status"] == "included"
        and feature["breakpoint_retained_percent"] == 100.0
        and feature["post_translation_status"] == "preserved"
        for feature in result["domains"]
    )


@pytest.mark.integration
def test_ewsr1_exon7_erg_exon6(reference_db: ReferenceDatabase) -> None:
    """
    EWSR1 exon 7 :: ERG exon 6 is a recurrent Ewing sarcoma fusion.

    The fixed GRCh38 breakpoints were chosen near the middle of the introns
    producing the exon 7::exon 6 junction:
      - 22:29287870 lies in EWSR1 intron 7/8 (22:29287135-29288605)
      - 21:38396509 lies in ERG intron 5/6 (21:38392445-38400573)

    This tests the same disruptive FET::ETS architecture with a different
    3' partner: the EWSR1 RNA-recognition motif (InterPro IPR000504) is
    excluded, while the ERG ETS domain (InterPro IPR000418) remains included
    and preserved.
    """
    result = _assert_success(
        annotate_fusion_domains(
            transcript1_id=GENE2TRANSCRIPT["EWSR1"],
            transcript2_id=GENE2TRANSCRIPT["ERG"],
            breakpoint1="22:29287870",
            breakpoint2="21:38396509",
            gene1_terminus="N",
            gene2_terminus="C",
            reference=reference_db,
        )
    )

    assert result["frame_status"] == "in_frame"

    assert any(
        feature["transcript_id"] == GENE2TRANSCRIPT["EWSR1"]
        and feature["interpro_id"] == "IPR000504"
        and feature["breakpoint_based_status"] == "excluded"
        and feature["breakpoint_retained_percent"] == 0.0
        for feature in result["domains"]
    )

    assert any(
        feature["transcript_id"] == GENE2TRANSCRIPT["ERG"]
        and feature["interpro_id"] == "IPR000418"
        and feature["breakpoint_based_status"] == "included"
        and feature["breakpoint_retained_percent"] == 100.0
        and feature["post_translation_status"] == "preserved"
        for feature in result["domains"]
    )
