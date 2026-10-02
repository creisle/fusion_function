from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal, NotRequired, TypedDict

from .reference import get_reference

if TYPE_CHECKING:
    from .data import ReferenceReader


Strand = Literal[-1, 1]


class ProteinFeatureAnnotation(TypedDict):
    name: str | None
    entry_type: str | None
    interpro_id: str | None


class TranscriptExon(TypedDict):
    exon_number: int
    chromosome: str | None
    genomic_start: int
    genomic_end: int
    premrna_start: int
    premrna_end: int


class CDSBlock(TypedDict):
    chromosome: str | None
    genomic_start: int
    genomic_end: int
    premrna_start: int
    premrna_end: int
    strand: Strand
    assembly_name: str | None
    cds_start: int
    cds_end: int


class GenomicSegment(TypedDict):
    chromosome: str | None
    start: int
    end: int
    strand: Strand
    assembly_name: str | None


class FeatureEvidence(TypedDict):
    code: str
    source: NotRequired[str]
    id: NotRequired[str]


class ProteinFeature(TypedDict):
    feature_id: str | None
    source: str | None
    interpro_id: str | None
    start: int
    end: int
    cds_start: int
    cds_end: int
    chromosome: str | None
    genomic_start: int | None
    genomic_end: int | None
    strand: Strand | None
    assembly_name: str | None
    description: str | None
    panther_subfamily_id: str | None
    panther_subfamily_description: str | None
    interpro_name: NotRequired[str | None]
    interpro_entry_type: NotRequired[str | None]
    genomic_segments: NotRequired[list[GenomicSegment]]
    feature_type: NotRequired[str]
    uniprot_accession: NotRequired[str]
    uniprot_isoform: NotRequired[str]
    evidence: NotRequired[list[FeatureEvidence]]


def feature_identity(feature: ProteinFeature) -> tuple[object, ...]:
    """Collapse only identical feature intervals, never overlapping boundaries.

    Integrated signatures may share an InterPro ID and exact interval. Features
    without one need their source, accession and description to distinguish, for
    example, different ligand-binding annotations at the same residue.
    """
    identity = (
        (feature["interpro_id"],)
        if feature["interpro_id"]
        else (
            feature["source"],
            feature["feature_id"],
            feature["description"],
            feature.get("uniprot_accession"),
            feature.get("uniprot_isoform"),
        )
    )
    return (
        *identity,
        feature["start"],
        feature["end"],
        feature.get("feature_type") or feature.get("interpro_entry_type"),
    )


class ReferenceSpliceSite(TypedDict):
    type: Literal["donor", "acceptor"]
    exon_number: int
    genomic_position: int
    premrna_position: int
    disruption_start: int
    disruption_end: int


class TranscriptProteinFeatureResult(TypedDict):
    translation_id: str | None
    protein_length: int | None
    chromosome: str | None
    strand: Strand
    transcript_genomic_start: int
    transcript_genomic_end: int
    premrna_length: int
    premrna_sequence: NotRequired[str]
    transcript_exons: list[TranscriptExon]
    cds_blocks: list[CDSBlock]
    cds_start_phase: NotRequired[int]
    protein_features: list[ProteinFeature]
    splice_sites: list[ReferenceSpliceSite]
    assembly_name: NotRequired[str | None]


class EnsemblError(TypedDict):
    error: str


ProteinFeatureResponse = TranscriptProteinFeatureResult | EnsemblError


def get_protein_domains(
    transcript_id: str,
    *,
    reference: ReferenceReader | None = None,
    database: str | Path | None = None,
    release: int | None = None,
) -> ProteinFeatureResponse:
    """Read prepared transcript structure, sequence and protein features locally."""
    reader = get_reference(reference=reference, database=database, release=release)
    return reader.get_transcript(transcript_id)
