from __future__ import annotations

from itertools import product
from typing import Literal, NotRequired, TypedDict, cast

from pathlib import Path

from .reference import ReferenceDatabase, get_reference

from .ensembl import (
    EnsemblError,
    ProteinFeature,
    Strand,
    TranscriptProteinFeatureResult,
    get_protein_domains,
)


FusionSideName = Literal["transcript1", "transcript2"]
FusionTerminus = Literal["N", "C"]
SpliceSiteType = Literal["donor", "acceptor"]
BreakpointBasedStatus = Literal["included", "disrupted", "excluded"]
PostSplicingStatus = Literal["preserved", "lost"]
PostTranslationStatus = Literal[
    "preserved",
    "frame_disrupted",
    "translation_start_disrupted",
    "translation_start_excluded",
    "premature_termination_disrupted",
    "premature_termination_excluded",
]
FrameStatus = Literal["in_frame", "out_of_frame"]
FusionFrameStatus = Literal["in_frame", "out_of_frame", "in_frame/out_of_frame"]


TranslationStartStatus = Literal[
    "native_start_retained",
    "alternative_start_found",
    "alternative_start_not_found",
    "native_start_lost",
]
TRANSLATION_START_STATUS_ORDER = (
    "native_start_retained",
    "alternative_start_found",
    "alternative_start_not_found",
    "native_start_lost",
)


class ProductSpliceSite(TypedDict):
    type: SpliceSiteType
    side: FusionSideName
    transcript_id: str
    exon_number: int


class CodingFragment(TypedDict):
    side: FusionSideName
    transcript_id: str
    source_cds_start: int
    source_cds_end: int
    product_start: int
    product_end: int


class SpliceProduct(TypedDict):
    product_index: int
    transcript_length: int
    translation_start: int | None
    stop_position: int | None
    frame_status: FrameStatus
    splice_sites: list[ProductSpliceSite]
    coding_fragments: list[CodingFragment]
    translation_start_source: NotRequired[Literal["alternative", "not_found"]]
    start_status: NotRequired[TranslationStartStatus]
    start_transcript_id: NotRequired[str | None]


class FusionCodingSegment(TypedDict):
    side: FusionSideName
    transcript_id: str
    fusion_terminus: FusionTerminus
    local_start: int
    local_end: int
    source_cds_start: int
    source_cds_end: int


class FusionSpliceSite(TypedDict):
    type: SpliceSiteType
    side: FusionSideName
    transcript_id: str
    exon_number: int
    genomic_position: int
    premrna_position: int
    local_position: int
    intact: bool


class FusionSide(TypedDict):
    side: FusionSideName
    transcript_id: str
    breakpoint: str
    fusion_terminus: FusionTerminus
    strand: Strand
    inferred_breakpoint_orientation: Literal["L", "R"]
    length: int
    sequence: str
    splice_sites: list[FusionSpliceSite]
    coding_segments: list[FusionCodingSegment]
    protein_features: list[ProteinFeature]


class LayoutSpliceSite(FusionSpliceSite):
    fusion_position: int


class LayoutCodingSegment(FusionCodingSegment):
    fusion_start: int
    fusion_end: int


class FusionLayout(TypedDict):
    n_terminal_side: FusionSideName
    total_length: int
    sequence: str
    splice_sites: list[LayoutSpliceSite]
    coding_segments: list[LayoutCodingSegment]


class FeatureProductPrediction(TypedDict):
    post_splicing_status: PostSplicingStatus
    post_translation_status: PostTranslationStatus | None


class AnnotatedProteinFeature(ProteinFeature):
    transcript_id: str
    breakpoint_based_status: BreakpointBasedStatus
    breakpoint_retained_percent: float
    post_splicing_status: str | None
    post_translation_status: str | None


class ResolvedAnnotatedProteinFeature(AnnotatedProteinFeature):
    resolved_interpro_id: str | None
    member_name: str | None
    member_entry_type: str | None


class AggregatedDomain(TypedDict):
    transcript_id: str
    interpro_id: str | None
    name: str | None
    specific_name: NotRequired[str | None]
    domain_type: str | None
    sources: list[str]
    start: int
    end: int
    breakpoint_based_status: str
    breakpoint_retained_percent: float
    post_splicing_status: str | None
    post_translation_status: str | None


class FunctionalDomain(TypedDict):
    transcript_id: str
    interpro_id: str | None
    name: str
    domain_type: str
    start: int
    end: int
    breakpoint_based_status: str
    breakpoint_retained_percent: float
    post_splicing_status: str | None
    post_translation_status: str | None


class FusionDomainResult(TypedDict):
    frame_status: FusionFrameStatus | None
    domains: list[FunctionalDomain]
    assumed_disruptive: NotRequired[bool]
    translation_start: NotRequired[dict[str, str]]


STOP_CODONS = {"TAA", "TAG", "TGA"}
RELEVANT_ENTRY_TYPES = {"domain", "binding_site", "active_site", "conserved_site"}


def _parse_breakpoint(breakpoint: str) -> tuple[str, int]:
    """Parse an exact chromosome:position breakpoint."""
    chromosome, position = breakpoint.split(":", 1)
    if "-" in position:
        raise NotImplementedError(f"Breakpoint intervals are not supported: {breakpoint}")
    return chromosome, int(position)


def _infer_breakpoint_orientation(
    strand: Strand, fusion_terminus: FusionTerminus
) -> Literal["L", "R"]:
    """Infer the retained genomic side from strand and fusion terminus."""
    return "L" if (strand, fusion_terminus) in {(1, "N"), (-1, "C")} else "R"


def _genomic_to_premrna(position: int, transcript: TranscriptProteinFeatureResult) -> int:
    """Convert a genomic position to transcript-oriented pre-mRNA coordinates."""
    start = transcript["transcript_genomic_start"]
    end = transcript["transcript_genomic_end"]
    if not start <= position <= end:
        raise ValueError(f"Position {position} lies outside transcript genomic span {start}-{end}")
    return position - start + 1 if transcript["strand"] == 1 else end - position + 1


def _build_fusion_side(
    side: FusionSideName,
    transcript_id: str,
    breakpoint: str,
    fusion_terminus: FusionTerminus,
    transcript: TranscriptProteinFeatureResult,
) -> FusionSide:
    """Build the retained pre-mRNA portion of one fusion partner."""
    chromosome, breakpoint_position = _parse_breakpoint(breakpoint)
    if str(transcript["chromosome"]).removeprefix("chr") != chromosome.removeprefix("chr"):
        raise ValueError(
            f"Breakpoint chromosome {chromosome} does not match "
            f"transcript chromosome {transcript['chromosome']}"
        )
    breakpoint_premrna = _genomic_to_premrna(breakpoint_position, transcript)
    if fusion_terminus == "N":
        retained_start, retained_end = 1, breakpoint_premrna
    else:
        retained_start, retained_end = breakpoint_premrna, transcript["premrna_length"]
    splice_sites: list[FusionSpliceSite] = []
    for site in transcript["splice_sites"]:
        premrna_position = site["premrna_position"]
        if retained_start <= premrna_position <= retained_end:
            splice_sites.append(
                {
                    "type": site["type"],
                    "side": side,
                    "transcript_id": transcript_id,
                    "exon_number": site["exon_number"],
                    "genomic_position": site["genomic_position"],
                    "premrna_position": premrna_position,
                    "local_position": premrna_position - retained_start + 1,
                    "intact": not site["disruption_start"]
                    <= breakpoint_position
                    <= site["disruption_end"],
                }
            )
    coding_segments: list[FusionCodingSegment] = []
    for block in transcript["cds_blocks"]:
        start = max(block["premrna_start"], retained_start)
        end = min(block["premrna_end"], retained_end)
        if start > end:
            continue
        coding_segments.append(
            {
                "side": side,
                "transcript_id": transcript_id,
                "fusion_terminus": fusion_terminus,
                "local_start": start - retained_start + 1,
                "local_end": end - retained_start + 1,
                "source_cds_start": block["cds_start"] + start - block["premrna_start"],
                "source_cds_end": block["cds_start"] + end - block["premrna_start"],
            }
        )
    return {
        "side": side,
        "transcript_id": transcript_id,
        "breakpoint": f"{chromosome}:{breakpoint_position}",
        "fusion_terminus": fusion_terminus,
        "strand": transcript["strand"],
        "inferred_breakpoint_orientation": _infer_breakpoint_orientation(
            transcript["strand"], fusion_terminus
        ),
        "length": retained_end - retained_start + 1,
        "sequence": transcript["premrna_sequence"][retained_start - 1 : retained_end],
        "splice_sites": splice_sites,
        "coding_segments": coding_segments,
        "protein_features": transcript["protein_features"],
    }


def _build_fusion_layout(
    side1: FusionSide, side2: FusionSide, inserted_sequence: str = ""
) -> FusionLayout:
    """Arrange retained partners N-to-C and build the fusion pre-mRNA."""
    if {side1["fusion_terminus"], side2["fusion_terminus"]} != {"N", "C"}:
        raise ValueError("Fusion prediction requires one N-terminal and one C-terminal partner")
    n_side, c_side = (side1, side2) if side1["fusion_terminus"] == "N" else (side2, side1)
    inserted_sequence = inserted_sequence.upper()
    c_offset = n_side["length"] + len(inserted_sequence)
    splice_sites: list[LayoutSpliceSite] = [
        {
            **site,
            "fusion_position": site["local_position"]
            + (0 if side["fusion_terminus"] == "N" else c_offset),
        }
        for side in (n_side, c_side)
        for site in side["splice_sites"]
    ]
    coding_segments: list[LayoutCodingSegment] = [
        {
            **segment,
            "fusion_start": segment["local_start"]
            + (0 if side["fusion_terminus"] == "N" else c_offset),
            "fusion_end": segment["local_end"]
            + (0 if side["fusion_terminus"] == "N" else c_offset),
        }
        for side in (n_side, c_side)
        for segment in side["coding_segments"]
    ]
    return {
        "n_terminal_side": n_side["side"],
        "total_length": n_side["length"] + len(inserted_sequence) + c_side["length"],
        "sequence": n_side["sequence"] + inserted_sequence + c_side["sequence"],
        "splice_sites": sorted(splice_sites, key=lambda x: x["fusion_position"]),
        "coding_segments": sorted(coding_segments, key=lambda x: x["fusion_start"]),
    }


def _build_retained_side_layout(side: FusionSide) -> FusionLayout:
    """Model the known retained sequence without inventing the unknown partner."""
    return {
        "n_terminal_side": side["side"],
        "total_length": side["length"],
        "sequence": side["sequence"],
        "splice_sites": [
            {**site, "fusion_position": site["local_position"]} for site in side["splice_sites"]
        ],
        "coding_segments": [
            {**segment, "fusion_start": segment["local_start"], "fusion_end": segment["local_end"]}
            for segment in side["coding_segments"]
        ],
    }


def _generate_splice_patterns(splice_sites: list[LayoutSpliceSite]) -> list[list[LayoutSpliceSite]]:
    """Generate splice patterns from intact donor/acceptor alternatives."""
    groups: list[list[LayoutSpliceSite]] = []
    for site in sorted(
        (site for site in splice_sites if site["intact"]), key=lambda x: x["fusion_position"]
    ):
        if groups and groups[-1][0]["type"] == site["type"]:
            groups[-1].append(site)
        else:
            groups.append([site])
    if groups and groups[0][0]["type"] == "acceptor":
        groups.pop(0)
    if groups and groups[-1][0]["type"] == "donor":
        groups.pop()
    return [list(pattern) for pattern in product(*groups)] if groups else [[]]


def _subtract_intervals(
    start: int, end: int, intervals: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Subtract inclusive intervals from another inclusive interval."""
    segments = [(start, end)]
    for remove_start, remove_end in intervals:
        updated = []
        for segment_start, segment_end in segments:
            if remove_end < segment_start or remove_start > segment_end:
                updated.append((segment_start, segment_end))
                continue
            if segment_start < remove_start:
                updated.append((segment_start, remove_start - 1))
            if segment_end > remove_end:
                updated.append((remove_end + 1, segment_end))
        segments = updated
    return segments


def _get_product_frame_status(
    translation_start: int | None, coding_fragments: list[CodingFragment]
) -> FrameStatus:
    """Determine whether retained CDS fragments preserve source codon phase."""
    if translation_start is None:
        return "out_of_frame"
    translated_fragments = [
        fragment for fragment in coding_fragments if fragment["product_end"] >= translation_start
    ]
    if not translated_fragments:
        return "out_of_frame"
    return (
        "in_frame"
        if all(
            (fragment["product_start"] - translation_start) % 3
            == (fragment["source_cds_start"] - 1) % 3
            for fragment in translated_fragments
        )
        else "out_of_frame"
    )


def _find_stop_position(sequence: str, translation_start: int | None) -> int | None:
    """Return the first base of the first in-frame stop codon."""
    if translation_start is None:
        return None
    for position in range(translation_start - 1, len(sequence) - 2, 3):
        if sequence[position : position + 3] in STOP_CODONS:
            return position + 1
    return None


def _build_splice_product(
    product_index: int,
    layout: FusionLayout,
    splice_pattern: list[LayoutSpliceSite],
    *,
    search_for_start: bool = False,
    validate_native_start: bool = False,
) -> SpliceProduct:
    """Construct one spliced fusion transcript and determine frame and stop."""
    if len(splice_pattern) % 2:
        raise ValueError("Splice pattern must contain donor/acceptor pairs")
    removed: list[tuple[int, int]] = []
    for donor, acceptor in zip(splice_pattern[::2], splice_pattern[1::2]):
        if donor["type"] != "donor" or acceptor["type"] != "acceptor":
            raise ValueError("Splice pattern must alternate donor and acceptor sites")
        if donor["fusion_position"] >= acceptor["fusion_position"]:
            raise ValueError("Donor must precede acceptor")
        if donor["fusion_position"] + 1 <= acceptor["fusion_position"] - 1:
            removed.append((donor["fusion_position"] + 1, acceptor["fusion_position"] - 1))
    retained: list[tuple[int, int, int, int]] = []
    product_position = 1
    for fusion_start, fusion_end in _subtract_intervals(1, layout["total_length"], removed):
        length = fusion_end - fusion_start + 1
        retained.append((fusion_start, fusion_end, product_position, product_position + length - 1))
        product_position += length
    sequence = "".join(
        layout["sequence"][fusion_start - 1 : fusion_end]
        for fusion_start, fusion_end, _, _ in retained
    )
    coding_fragments: list[CodingFragment] = []
    for segment in layout["coding_segments"]:
        for fusion_start, fusion_end, product_start, _ in retained:
            start = max(segment["fusion_start"], fusion_start)
            end = min(segment["fusion_end"], fusion_end)
            if start > end:
                continue
            coding_fragments.append(
                {
                    "side": segment["side"],
                    "transcript_id": segment["transcript_id"],
                    "source_cds_start": segment["source_cds_start"]
                    + start
                    - segment["fusion_start"],
                    "source_cds_end": segment["source_cds_start"] + end - segment["fusion_start"],
                    "product_start": product_start + start - fusion_start,
                    "product_end": product_start + end - fusion_start,
                }
            )
    coding_fragments.sort(key=lambda x: x["product_start"])
    # Select one native initiation site in the complete spliced product. A
    # retained C-terminal start can supply initiation if the N portion is UTR-only.
    native_candidates = [
        (fragment["product_start"] + 1 - fragment["source_cds_start"], fragment["transcript_id"])
        for fragment in coding_fragments
        if fragment["source_cds_start"] <= 1 <= fragment["source_cds_end"]
    ]
    if search_for_start or validate_native_start:
        native_candidates = [
            (position, transcript_id)
            for position, transcript_id in native_candidates
            if sequence[position - 1 : position + 2] == "ATG"
        ]
    translation_start, start_transcript_id = next(iter(native_candidates), (None, None))
    start_source: Literal["alternative", "not_found"] | None = None
    native_start_retained = translation_start is not None
    if translation_start is None and search_for_start and coding_fragments:
        # Search after splicing, in any phase. Choosing the earliest retained ATG
        # is a heuristic; compatibility with the original protein is checked
        # separately. With no retained CDS there are no original domains to translate.
        start_index = sequence.find("ATG")
        translation_start = start_index + 1 if start_index >= 0 else None
        start_source = "alternative" if translation_start is not None else "not_found"
    result: SpliceProduct = {
        "product_index": product_index,
        "transcript_length": product_position - 1,
        "translation_start": translation_start,
        "stop_position": _find_stop_position(sequence, translation_start),
        "frame_status": _get_product_frame_status(translation_start, coding_fragments),
        "splice_sites": [
            {
                "type": site["type"],
                "side": site["side"],
                "transcript_id": site["transcript_id"],
                "exon_number": site["exon_number"],
            }
            for site in splice_pattern
        ],
        "coding_fragments": coding_fragments,
    }
    if start_source is not None:
        result["translation_start_source"] = start_source
    # This describes initiation of this product once, rather than separate
    # outcomes for each occurrence of the same source transcript.
    result["start_status"] = (
        "native_start_retained"
        if native_start_retained
        else "native_start_lost"
        if not coding_fragments or not search_for_start
        else "alternative_start_found"
        if translation_start is not None
        else "alternative_start_not_found"
    )
    result["start_transcript_id"] = start_transcript_id
    return result


def _covered_length(start: int, end: int, intervals: list[tuple[int, int]]) -> int:
    """Return the number of bases covered within an inclusive interval."""
    intervals = sorted(
        (max(start, a), min(end, b)) for a, b in intervals if b >= start and a <= end
    )
    if not intervals:
        return 0
    covered = 0
    current_start, current_end = intervals[0]
    for interval_start, interval_end in intervals[1:]:
        if interval_start <= current_end + 1:
            current_end = max(current_end, interval_end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = interval_start, interval_end
    return covered + current_end - current_start + 1


def _classify_feature(
    feature: ProteinFeature, side: FusionSideName, splice_product: SpliceProduct
) -> FeatureProductPrediction:
    """Check the feature portion retained by the breakpoint in one splice product.

    The caller clips CDS bounds to the retained portion. This lets a disrupted
    domain undergo the same downstream checks as a fully included domain.
    """
    fragments: list[tuple[int, int, int, int]] = []
    for fragment in splice_product["coding_fragments"]:
        if (
            fragment["side"] != side
            or fragment["source_cds_end"] < feature["cds_start"]
            or fragment["source_cds_start"] > feature["cds_end"]
        ):
            continue
        source_start = max(feature["cds_start"], fragment["source_cds_start"])
        source_end = min(feature["cds_end"], fragment["source_cds_end"])
        fragments.append(
            (
                source_start,
                source_end,
                fragment["product_start"] + source_start - fragment["source_cds_start"],
                fragment["product_start"] + source_end - fragment["source_cds_start"],
            )
        )
    fragments.sort()
    if (
        _covered_length(
            feature["cds_start"],
            feature["cds_end"],
            [(start, end) for start, end, _, _ in fragments],
        )
        < feature["cds_end"] - feature["cds_start"] + 1
    ):
        return {"post_splicing_status": "lost", "post_translation_status": None}
    translation_start = splice_product["translation_start"]
    if translation_start is None:
        status: PostTranslationStatus = (
            "translation_start_excluded"
            if splice_product.get("translation_start_source") == "not_found"
            else "frame_disrupted"
        )
    elif translation_start > fragments[-1][3]:
        status = "translation_start_excluded"
    else:
        # An alternative ATG can occur inside a retained feature. Compare phase
        # and continuity only for bases translated from that start onward.
        translated = [
            (
                start + max(0, translation_start - product_start),
                end,
                max(product_start, translation_start),
                product_end,
            )
            for start, end, product_start, product_end in fragments
            if product_end >= translation_start
        ]
        if (
            any(
                current[0] != previous[1] + 1 or current[2] != previous[3] + 1
                for previous, current in zip(translated, translated[1:])
            )
            or (translated[0][2] - translation_start) % 3 != (translated[0][0] - 1) % 3
        ):
            status = "frame_disrupted"
        elif (
            splice_product["stop_position"] is not None
            and splice_product["stop_position"] <= translated[0][2]
        ):
            status = "premature_termination_excluded"
        elif (
            splice_product["stop_position"] is not None
            and splice_product["stop_position"] <= translated[-1][3]
        ):
            status = "premature_termination_disrupted"
        elif translation_start > fragments[0][2]:
            status = "translation_start_disrupted"
        else:
            status = "preserved"
    return {"post_splicing_status": "preserved", "post_translation_status": status}


def _combine_statuses(statuses: list[str | None], order: tuple[str, ...]) -> str | None:
    """Combine statuses into a stable slash-delimited summary."""
    observed = {status for value in statuses if value is not None for status in value.split("/")}
    return "/".join(status for status in order if status in observed) or None


def _annotate_feature_statuses(
    side: FusionSide, splice_products: list[SpliceProduct]
) -> list[AnnotatedProteinFeature]:
    """Annotate each protein feature independently across predicted products."""
    retained_cds = [
        (segment["source_cds_start"], segment["source_cds_end"])
        for segment in side["coding_segments"]
    ]
    annotated: list[AnnotatedProteinFeature] = []
    for feature in side["protein_features"]:
        feature_length = feature["cds_end"] - feature["cds_start"] + 1
        covered = _covered_length(feature["cds_start"], feature["cds_end"], retained_cds)
        breakpoint_status: BreakpointBasedStatus = (
            "excluded" if covered == 0 else "disrupted" if covered < feature_length else "included"
        )
        if covered > 0:
            # An N- or C-terminal breakpoint retains a continuous source CDS
            # interval, even when its genomic sequence spans several exons.
            # Check that interval after splicing without treating bases already
            # removed by the breakpoint as an additional splicing loss.
            retained_feature: ProteinFeature = {
                **feature,
                "cds_start": max(feature["cds_start"], min(start for start, _ in retained_cds)),
                "cds_end": min(feature["cds_end"], max(end for _, end in retained_cds)),
            }
            predictions = [
                _classify_feature(retained_feature, side["side"], splice_product)
                for splice_product in splice_products
            ]
            splicing_status = _combine_statuses(
                [x["post_splicing_status"] for x in predictions], ("preserved", "lost")
            )
            sequence_status = _combine_statuses(
                [x["post_translation_status"] for x in predictions],
                (
                    "preserved",
                    "frame_disrupted",
                    "translation_start_disrupted",
                    "translation_start_excluded",
                    "premature_termination_disrupted",
                    "premature_termination_excluded",
                ),
            )
        else:
            splicing_status = None
            sequence_status = None
        annotated.append(
            {
                **feature,
                "transcript_id": side["transcript_id"],
                "breakpoint_based_status": breakpoint_status,
                "breakpoint_retained_percent": round(100 * covered / feature_length, 1),
                "post_splicing_status": splicing_status,
                "post_translation_status": sequence_status,
            }
        )
    return annotated


def _resolve_feature_metadata(
    features: list[AnnotatedProteinFeature],
) -> list[ResolvedAnnotatedProteinFeature]:
    """Use metadata already resolved during preprocessing."""
    resolved: list[ResolvedAnnotatedProteinFeature] = []
    for feature in features:
        entry_type = feature.get("interpro_entry_type")
        if feature["interpro_id"] and entry_type is None:
            raise ValueError(
                f"Missing InterPro entry type for {feature['interpro_id']}. "
                "Reprocess the database with 'fusion-function prepare-data --preprocess-only DB'. "
                "Metadata is fetched automatically; --interpro-entries FILE supplies a local list."
            )
        resolved.append(
            {
                **feature,
                "resolved_interpro_id": feature["interpro_id"],
                "member_name": feature.get("interpro_name"),
                "member_entry_type": entry_type,
            }
        )
    return resolved


def _overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """Return overlap between two inclusive intervals."""
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def _aggregate_feature_group(features: list[ResolvedAnnotatedProteinFeature]) -> AggregatedDomain:
    """Aggregate feature records representing one InterPro occurrence."""
    interpro_id = next(
        (
            feature["resolved_interpro_id"]
            for feature in features
            if feature["resolved_interpro_id"] is not None
        ),
        None,
    )
    return {
        "transcript_id": features[0]["transcript_id"],
        "interpro_id": interpro_id,
        "name": next(
            (feature["member_name"] for feature in features if feature["member_name"]),
            next((feature["description"] for feature in features if feature["description"]), None),
        ),
        "specific_name": next(
            (
                feature.get("panther_subfamily_description")
                for feature in features
                if feature.get("panther_subfamily_description")
            ),
            None,
        ),
        "domain_type": next(
            (feature["member_entry_type"] for feature in features if feature["member_entry_type"]),
            None,
        ),
        "sources": sorted(
            {feature["source"] for feature in features if feature["source"] is not None}
        ),
        "start": min(feature["start"] for feature in features),
        "end": max(feature["end"] for feature in features),
        "breakpoint_based_status": cast(
            str,
            _combine_statuses(
                [feature["breakpoint_based_status"] for feature in features],
                ("included", "disrupted", "excluded"),
            ),
        ),
        "breakpoint_retained_percent": max(
            feature["breakpoint_retained_percent"] for feature in features
        ),
        "post_splicing_status": _combine_statuses(
            [feature["post_splicing_status"] for feature in features], ("preserved", "lost")
        ),
        "post_translation_status": _combine_statuses(
            [feature["post_translation_status"] for feature in features],
            (
                "preserved",
                "frame_disrupted",
                "translation_start_disrupted",
                "translation_start_excluded",
                "premature_termination_disrupted",
                "premature_termination_excluded",
            ),
        ),
    }


def _same_functional_feature(feature1: AggregatedDomain, feature2: AggregatedDomain) -> bool:
    """Determine whether two annotations represent one functional feature."""
    if (
        feature1["transcript_id"] != feature2["transcript_id"]
        or feature1["domain_type"] != feature2["domain_type"]
    ):
        return False
    overlap = _overlap(feature1["start"], feature1["end"], feature2["start"], feature2["end"])
    if overlap == 0:
        return False
    length1 = feature1["end"] - feature1["start"] + 1
    length2 = feature2["end"] - feature2["start"] + 1
    if (feature1["interpro_id"] is None) != (feature2["interpro_id"] is None):
        unresolved_length = length1 if feature1["interpro_id"] is None else length2
        return overlap / unresolved_length >= 0.7
    return min(overlap / length1, overlap / length2) >= 0.7


def _finalize_domains(domains: list[AggregatedDomain]) -> list[FunctionalDomain]:
    """Filter and collapse redundant functional protein annotations."""

    # Non-functional annotations with a useful protein/family-specific name.
    # These are retained only as fallbacks when no real functional InterPro
    # feature represents the same region.
    specific_annotations = [
        domain
        for domain in domains
        if domain.get("specific_name") and domain["domain_type"] not in RELEVANT_ENTRY_TYPES
    ]

    retained_domains = [
        cast(AggregatedDomain, dict(domain))
        for domain in domains
        if domain["domain_type"] in RELEVANT_ENTRY_TYPES
    ]

    # Suppress specific family annotations that are already represented by
    # a real InterPro domain. Do NOT copy specific_name onto the domain:
    # the canonical InterPro name should remain the functional feature name.
    represented_specific_annotations: set[int] = set()

    for domain in retained_domains:
        if domain["domain_type"] != "domain":
            continue

        domain_length = domain["end"] - domain["start"] + 1

        for index, annotation in enumerate(specific_annotations):
            if annotation["transcript_id"] != domain["transcript_id"]:
                continue

            overlap = _overlap(
                domain["start"], domain["end"], annotation["start"], annotation["end"]
            )
            annotation_length = annotation["end"] - annotation["start"] + 1

            if overlap / min(domain_length, annotation_length) >= 0.7:
                represented_specific_annotations.add(index)

    # If a specific family annotation is not represented by an actual
    # functional feature, retain it as a fallback domain. In this case the
    # specific_name is useful because there is no canonical domain name.
    for index, annotation in enumerate(specific_annotations):
        if index in represented_specific_annotations:
            continue

        retained_domains.append(
            {**annotation, "name": cast(str, annotation["specific_name"]), "domain_type": "domain"}
        )

    buckets: dict[tuple[str, str], list[AggregatedDomain]] = {}

    for domain in retained_domains:
        buckets.setdefault((domain["transcript_id"], cast(str, domain["domain_type"])), []).append(
            domain
        )

    final: list[FunctionalDomain] = []

    for bucket in buckets.values():
        groups: list[list[AggregatedDomain]] = []

        for domain in sorted(
            bucket,
            key=lambda x: (x["interpro_id"] is not None, len(x["sources"]), x["end"] - x["start"]),
            reverse=True,
        ):
            for group in groups:
                if _same_functional_feature(domain, group[0]):
                    group.append(domain)
                    break
            else:
                groups.append([domain])

        for group in groups:
            representative = group[0]

            final.append(
                {
                    "transcript_id": representative["transcript_id"],
                    "interpro_id": representative["interpro_id"],
                    "name": (
                        representative["name"]
                        or representative["interpro_id"]
                        or "Unnamed functional feature"
                    ),
                    "domain_type": cast(str, representative["domain_type"]),
                    "start": representative["start"],
                    "end": representative["end"],
                    "breakpoint_based_status": cast(
                        str,
                        _combine_statuses(
                            [x["breakpoint_based_status"] for x in group],
                            ("included", "disrupted", "excluded"),
                        ),
                    ),
                    "breakpoint_retained_percent": representative["breakpoint_retained_percent"],
                    "post_splicing_status": _combine_statuses(
                        [x["post_splicing_status"] for x in group], ("preserved", "lost")
                    ),
                    "post_translation_status": _combine_statuses(
                        [x["post_translation_status"] for x in group],
                        (
                            "preserved",
                            "frame_disrupted",
                            "translation_start_disrupted",
                            "translation_start_excluded",
                            "premature_termination_disrupted",
                            "premature_termination_excluded",
                        ),
                    ),
                }
            )

    return sorted(
        final, key=lambda x: (x["transcript_id"], x["start"], x["end"], x["domain_type"], x["name"])
    )


def _summarize_frame_status(splice_products: list[SpliceProduct]) -> FusionFrameStatus:
    """Summarize frame status across predicted splice products."""
    statuses = {product["frame_status"] for product in splice_products}
    if statuses == {"in_frame"}:
        return "in_frame"
    if statuses == {"out_of_frame"}:
        return "out_of_frame"
    return "in_frame/out_of_frame"


def _summarize_translation_start(products: list[SpliceProduct]) -> str:
    """Combine distinct product initiation outcomes in a stable slash order."""
    return cast(
        str,
        _combine_statuses(
            [product["start_status"] for product in products], TRANSLATION_START_STATUS_ORDER
        ),
    )


def annotate_fusion_domains(
    transcript1_id: str | None = None,
    transcript2_id: str | None = None,
    breakpoint1: str | None = None,
    breakpoint2: str | None = None,
    gene1_terminus: FusionTerminus | None = None,
    gene2_terminus: FusionTerminus | None = None,
    inserted_sequence: str = "",
    *,
    reference: ReferenceDatabase | None = None,
    database: str | Path | None = None,
    release: int | None = None,
) -> FusionDomainResult | EnsemblError:
    """Annotate retained features and reconstruct a fusion when both ends permit it.

    A single known partner is assumed disruptive and its retained portion is
    spliced and checked. Every product lacking a native start but retaining CDS
    uses its first ATG as a heuristic start. N/N and C/C pairs are unsupported.
    """
    partners = [
        ("transcript1", transcript1_id, breakpoint1, gene1_terminus),
        ("transcript2", transcript2_id, breakpoint2, gene2_terminus),
    ]
    if (
        transcript1_id
        and transcript2_id
        and gene1_terminus == gene2_terminus
        and gene1_terminus in ("N", "C")
    ):
        raise NotImplementedError(
            "N/N and C/C fusions are unsupported; two known partners require one N and one C terminus"
        )
    known_partners: list[tuple[FusionSideName, str, str, FusionTerminus]] = []
    for side_name, transcript_id, breakpoint, terminus in partners:
        if not transcript_id:
            continue
        if not breakpoint or terminus not in ("N", "C"):
            return {
                "error": f"{side_name} ({transcript_id}) requires a breakpoint and "
                f"an N or C terminus; received breakpoint={breakpoint!r}, terminus={terminus!r}"
            }
        known_partners.append(
            (
                cast(FusionSideName, side_name),
                transcript_id,
                breakpoint,
                cast(FusionTerminus, terminus),
            )
        )
    if not known_partners:
        return {"error": "At least one known transcript is required for feature annotation"}

    reader = get_reference(reference=reference, database=database, release=release)
    transcripts: dict[str, TranscriptProteinFeatureResult] = {}
    for transcript_id in {partner[1] for partner in known_partners}:
        transcript_result = get_protein_domains(transcript_id, reference=reader)
        if "error" in transcript_result:
            return cast(EnsemblError, transcript_result)
        transcripts[transcript_id] = cast(TranscriptProteinFeatureResult, transcript_result)
    sides: list[FusionSide] = []
    for side_name, transcript_id, breakpoint, terminus in known_partners:
        try:
            sides.append(
                _build_fusion_side(
                    side_name, transcript_id, breakpoint, terminus, transcripts[transcript_id]
                )
            )
        except (ValueError, NotImplementedError) as error:
            return {"error": f"{side_name} ({transcript_id}): {error}"}

    can_reconstruct = len(sides) == 2 and {side["fusion_terminus"] for side in sides} == {"N", "C"}
    splice_products: list[SpliceProduct] = []
    products_by_side: dict[FusionSideName, list[SpliceProduct]] = {}
    if can_reconstruct:
        layout = _build_fusion_layout(sides[0], sides[1], inserted_sequence)
        splice_products = [
            _build_splice_product(
                index, layout, pattern, search_for_start=True, validate_native_start=True
            )
            for index, pattern in enumerate(_generate_splice_patterns(layout["splice_sites"]), 1)
        ]
        products_by_side = {side["side"]: splice_products for side in sides}
    else:
        for side in sides:
            layout = _build_retained_side_layout(side)
            products = [
                _build_splice_product(
                    index, layout, pattern, search_for_start=True, validate_native_start=True
                )
                for index, pattern in enumerate(
                    _generate_splice_patterns(layout["splice_sites"]), 1
                )
            ]
            products_by_side[side["side"]] = products
            splice_products.extend(products)
    # Group membership is reference-only and has already been computed.
    # Merge matching groups when the same transcript is used on both sides.
    grouped: dict[tuple[str, int], list[ResolvedAnnotatedProteinFeature]] = {}
    for side in sides:
        transcript_id = side["transcript_id"]
        features = _resolve_feature_metadata(
            _annotate_feature_statuses(side, products_by_side[side["side"]])
        )
        for group_index, indices in enumerate(transcripts[transcript_id]["feature_groups"]):
            grouped.setdefault((transcript_id, group_index), []).extend(
                features[index] for index in indices
            )
    domains = [_aggregate_feature_group(features) for features in grouped.values()]
    result: FusionDomainResult = {
        "frame_status": (
            _summarize_frame_status(splice_products)
            if can_reconstruct
            or any(
                product["translation_start"] is not None
                and any(
                    fragment["product_end"] >= product["translation_start"]
                    for fragment in product["coding_fragments"]
                )
                for product in splice_products
            )
            else None
        ),
        "domains": _finalize_domains(domains),
    }
    # Each complete product has one initiation prediction. Key it by the
    # transcript supplying the selected native start, or the intended initiating
    # transcript when that start is lost or an alternative start is used.
    default_start_transcript = next(
        (side["transcript_id"] for side in sides if side["fusion_terminus"] == "N"),
        sides[0]["transcript_id"],
    )
    start_products: dict[str, list[SpliceProduct]] = {}
    for product in splice_products:
        transcript_id = product["start_transcript_id"] or default_start_transcript
        start_products.setdefault(transcript_id, []).append(product)
    result["translation_start"] = {
        transcript_id: _summarize_translation_start(products)
        for transcript_id, products in start_products.items()
    }
    if not can_reconstruct:
        result["assumed_disruptive"] = True
    return result
