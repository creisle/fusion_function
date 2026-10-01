from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal, NotRequired, TypedDict

from .reference import get_reference

if TYPE_CHECKING:
    from .data import ReferenceReader


Strand = Literal[-1, 1]


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
    protein_features: list[ProteinFeature]
    splice_sites: list[ReferenceSpliceSite]
    feature_groups: list[list[int]]
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
