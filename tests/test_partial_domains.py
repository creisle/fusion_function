"""Downstream checks must evaluate the surviving portion of a cut domain."""

from typing import cast

import pytest

from fusion_function.ensembl import ProteinFeature
from fusion_function.fusion import FusionSide, SpliceProduct, _annotate_feature_statuses


def partial_side(terminus="C"):
    """A ten-residue domain cut within a codon, spanning two CDS blocks."""
    intervals = [(1, 6), (7, 20)] if terminus == "N" else [(5, 12), (13, 30)]
    feature = cast(
        ProteinFeature,
        {
            "feature_id": "test-domain",
            "interpro_id": "IPR000001",
            "start": 1,
            "end": 10,
            "cds_start": 1,
            "cds_end": 30,
        },
    )
    return cast(
        FusionSide,
        {
            "side": "transcript1" if terminus == "N" else "transcript2",
            "transcript_id": "test-transcript",
            "fusion_terminus": terminus,
            "coding_segments": [
                {"source_cds_start": start, "source_cds_end": end} for start, end in intervals
            ],
            "protein_features": [feature],
        },
    )


def product_for(side, *, shift=0, stop=None, translation_start=1):
    """Map retained source CDS bases contiguously into a product in native phase."""
    first_source = side["coding_segments"][0]["source_cds_start"]
    product_start = 1 + (first_source - 1) % 3 + shift
    fragments = []
    for segment in side["coding_segments"]:
        start, end = segment["source_cds_start"], segment["source_cds_end"]
        fragments.append(
            {
                "side": side["side"],
                "transcript_id": side["transcript_id"],
                "source_cds_start": start,
                "source_cds_end": end,
                "product_start": product_start,
                "product_end": product_start + end - start,
            }
        )
        product_start += end - start + 1
    return cast(
        SpliceProduct,
        {
            "coding_fragments": fragments,
            "translation_start": translation_start,
            "stop_position": stop,
        },
    )


@pytest.mark.parametrize("terminus,percent", [("N", 66.7), ("C", 86.7)])
def test_partial_domains_preserve_retained_portion(terminus, percent):
    side = partial_side(terminus)
    original_feature = dict(side["protein_features"][0])
    (result,) = _annotate_feature_statuses(side, [product_for(side)])
    assert result["breakpoint_based_status"] == "disrupted"
    assert result["breakpoint_retained_percent"] == percent
    assert result["splicing_based_status"] == "preserved"
    assert result["fusion_sequence_status"] == "preserved"
    # Public bounds and the cached reference describe the original whole domain.
    assert (result["start"], result["end"]) == (1, 10)
    assert (result["cds_start"], result["cds_end"]) == (1, 30)
    assert side["protein_features"][0] == original_feature


@pytest.mark.parametrize("terminus", ["N", "C"])
@pytest.mark.parametrize(
    "stop,expected",
    [
        (None, "preserved"),
        (31, "preserved"),
        (1, "premature_termination_excluded"),
        (7, "premature_termination_disrupted"),
    ],
)
def test_partial_domains_check_premature_stops(terminus, stop, expected):
    side = partial_side(terminus)
    (result,) = _annotate_feature_statuses(side, [product_for(side, stop=stop)])
    assert result["splicing_based_status"] == "preserved"
    assert result["fusion_sequence_status"] == expected


@pytest.mark.parametrize("terminus", ["N", "C"])
@pytest.mark.parametrize("options", [{"shift": 1}, {"translation_start": None}])
def test_partial_domains_check_frame(terminus, options):
    side = partial_side(terminus)
    (result,) = _annotate_feature_statuses(side, [product_for(side, **options)])
    assert result["splicing_based_status"] == "preserved"
    assert result["fusion_sequence_status"] == "frame_disrupted"


@pytest.mark.parametrize("terminus", ["N", "C"])
@pytest.mark.parametrize("all_removed", [False, True])
def test_partial_domains_detect_further_splicing_loss(terminus, all_removed):
    side = partial_side(terminus)
    product = product_for(side)
    if all_removed:
        product["coding_fragments"] = []
    else:
        # Splicing skips the last retained exon, removing additional domain bases.
        product["coding_fragments"].pop()
    (result,) = _annotate_feature_statuses(side, [product])
    assert result["breakpoint_based_status"] == "disrupted"
    assert result["splicing_based_status"] == "lost"
    assert result["fusion_sequence_status"] is None


def test_partial_domain_combines_alternative_products():
    side = partial_side()
    lost = product_for(side)
    lost["coding_fragments"] = []
    (result,) = _annotate_feature_statuses(
        side, [product_for(side), product_for(side, shift=1), lost]
    )
    assert result["splicing_based_status"] == "preserved/lost"
    assert result["fusion_sequence_status"] == "preserved/frame_disrupted"


def test_excluded_domain_has_no_downstream_status():
    side = partial_side()
    side["coding_segments"] = [{"source_cds_start": 31, "source_cds_end": 60}]
    (result,) = _annotate_feature_statuses(side, [product_for(side)])
    assert result["breakpoint_based_status"] == "excluded"
    assert result["breakpoint_retained_percent"] == 0.0
    assert result["splicing_based_status"] is None
    assert result["fusion_sequence_status"] is None
