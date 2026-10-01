"""Local InterPro metadata access (no network requests)."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

from .reference import get_reference

if TYPE_CHECKING:
    from .data import ReferenceReader


class ProteinFeatureAnnotation(TypedDict):
    name: str | None
    entry_type: str | None
    interpro_id: str | None


def get_interpro_annotation(
    interpro_id: str, *, reference: ReferenceReader | None = None,
    database: str | Path | None = None, release: int | None = None,
) -> ProteinFeatureAnnotation | None:
    return get_reference(reference=reference, database=database, release=release).get_interpro_annotation(interpro_id)
