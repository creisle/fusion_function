"""Reviewed human UniProt features, imported once during reference preparation.

XML is streamed with the standard library. Isoforms are reconstructed only from
explicit UniProt splice variants; canonical features are transferred only across
unchanged sequence. Ensembl mappings additionally require an exact whole-protein
sequence match. No gene-name matching, approximate alignment or runtime network
access is used.
"""

from __future__ import annotations

import json
import sqlite3
from collections import OrderedDict
from collections.abc import Iterator
from pathlib import Path
from typing import TypedDict
from xml.etree import ElementTree as ET

from .ensembl import FeatureEvidence


UNIPROT_HUMAN_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=xml&compressed=true"
    "&query=organism_id%3A9606%20AND%20reviewed%3Atrue"
)
NS = "{http://uniprot.org/uniprot}"
FEATURE_TYPES = {
    "domain": "domain",
    "active site": "active_site",
    "binding site": "binding_site",
    "site": "conserved_site",
    "motif": "motif",
    "short sequence motif": "motif",
}


class CuratedFeature(TypedDict):
    feature_id: str
    feature_type: str
    start: int
    end: int
    description: str
    uniprot_accession: str
    uniprot_isoform: str
    evidence: list[FeatureEvidence]


class CuratedProtein(TypedDict):
    sequence: str
    features: list[CuratedFeature]
    accessions: list[str]
    protein_ids: list[str]


def exact_interval(element: ET.Element) -> tuple[int, int] | None:
    """Reject uncertain/fuzzy coordinates rather than inventing exact bounds."""
    location = element.find(NS + "location")
    if location is None:
        return None
    position = location.find(NS + "position")
    nodes = (
        [position, position]
        if position is not None
        else [location.find(NS + "begin"), location.find(NS + "end")]
    )
    if any(
        node is None
        or node.get("status") not in {None, "certain"}
        or not (node.get("position") or "").isdigit()
        for node in nodes
    ):
        return None
    start, end = [int(node.attrib["position"]) for node in nodes if node is not None]
    return (start, end) if 1 <= start <= end else None


def entry_proteins(entry: ET.Element) -> Iterator[CuratedProtein]:
    """Yield canonical/described isoforms with safely mapped functional features."""
    if entry.get("dataset") != "Swiss-Prot" or not any(
        ref.get("type") == "NCBI Taxonomy" and ref.get("id") == "9606"
        for ref in entry.findall(f"{NS}organism/{NS}dbReference")
    ):
        return
    accessions = [node.text for node in entry.findall(NS + "accession") if node.text]
    canonical = "".join((entry.findtext(NS + "sequence") or "").split())
    if not accessions or not canonical:
        raise ValueError("Reviewed human UniProt entry lacks accession or sequence")
    accession = accessions[0]
    evidence: dict[str, FeatureEvidence] = {}
    for node in entry.findall(NS + "evidence"):
        item: FeatureEvidence = {"code": node.get("type", "")}
        reference = node.find(f"{NS}source/{NS}dbReference")
        if reference is not None:
            item["source"] = reference.get("type", "")
            item["id"] = reference.get("id", "")
        evidence[node.get("key", "")] = item
    features = entry.findall(NS + "feature")
    variants = {node.get("id"): node for node in features if node.get("type") == "splice variant"}
    isoforms = entry.findall(f"{NS}comment[@type='alternative products']/{NS}isoform")
    # A canonical entry need not declare alternative products.
    sequences: dict[str, tuple[str, list[tuple[int, int, int]]]] = {accession: (canonical, [])}
    displayed = accession
    for isoform in isoforms:
        sequence_element = isoform.find(NS + "sequence")
        ids = [node.text for node in isoform.findall(NS + "id") if node.text]
        if sequence_element is None or not ids:
            continue
        edits: list[tuple[int, int, str]] = []
        if sequence_element.get("type") == "displayed":
            displayed = ids[0]
        elif sequence_element.get("type") == "described":
            refs = (sequence_element.get("ref") or "").split()
            if not refs:
                continue
            for variant_id in refs:
                variant = variants.get(variant_id)
                interval = exact_interval(variant) if variant is not None else None
                if variant is None or interval is None or interval[1] > len(canonical):
                    break
                start, end = interval
                original = variant.findtext(NS + "original")
                if original and original != canonical[start - 1 : end]:
                    break
                replacement = variant.findtext(NS + "variation") or ""
                edits.append((start, end, replacement))
            if len(edits) != len(refs):
                continue
        else:
            # External/not-described isoform sequences cannot be reconstructed.
            continue
        edits.sort()
        if any(left[1] >= right[0] for left, right in zip(edits, edits[1:])):
            continue
        parts: list[str] = []
        cursor = 0
        for start, end, replacement in edits:
            parts.extend((canonical[cursor : start - 1], replacement))
            cursor = end
        parts.append(canonical[cursor:])
        reconstructed = "".join(parts)
        if not reconstructed:
            continue
        for isoform_id in ids:
            sequences[isoform_id] = (reconstructed, [(a, b, len(c)) for a, b, c in edits])
    if displayed != accession:
        sequences.pop(accession)
    for isoform_id, (sequence, edits_map) in sequences.items():
        # The accession and displayed isoform are aliases for the same product.
        aliases = {isoform_id, accession} if isoform_id in {accession, displayed} else {isoform_id}
        protein_ids: list[str] = []
        for ref in entry.findall(f"{NS}dbReference[@type='Ensembl']"):
            molecule = ref.find(NS + "molecule")
            if (molecule.get("id") if molecule is not None else accession) not in aliases:
                continue
            protein_ids.extend(
                node.get("value", "").split(".")[0]
                for node in ref.findall(f"{NS}property[@type='protein sequence ID']")
            )
        curated: list[CuratedFeature] = []
        for index, node in enumerate(features, 1):
            feature_type = FEATURE_TYPES.get(node.get("type", ""))
            interval = exact_interval(node)
            if feature_type is None or interval is None:
                continue
            start, end = interval
            location = node.find(NS + "location")
            target = location.get("sequence") if location is not None else None
            if target:
                if target not in aliases:
                    continue
            else:
                # A feature touching edited residues is omitted, even when the
                # replacement has equal length. Unchanged features can shift.
                if any(start <= b and end >= a for a, b, _ in edits_map):
                    continue
                shift = sum(size - (b - a + 1) for a, b, size in edits_map if b < start)
                start, end = start + shift, end + shift
            if not 1 <= start <= end <= len(sequence):
                continue
            description = (
                node.get("description")
                or node.findtext(f"{NS}ligand/{NS}name")
                or node.get("type", "")
                or "Unnamed UniProt feature"
            )
            curated.append(
                {
                    "feature_id": node.get("id") or f"{accession}:{index}",
                    "feature_type": feature_type,
                    "start": start,
                    "end": end,
                    "description": description,
                    "uniprot_accession": accession,
                    "uniprot_isoform": displayed if isoform_id == accession else isoform_id,
                    "evidence": [
                        evidence[key]
                        for key in (node.get("evidence") or "").split()
                        if key in evidence
                    ],
                }
            )
        yield {
            "sequence": sequence,
            "features": curated,
            "accessions": list(
                aliases | (set(accessions) if isoform_id in {accession, displayed} else set())
            ),
            "protein_ids": protein_ids,
        }


def import_features(db: sqlite3.Connection, path: Path) -> dict[str, int]:
    """Store sequence-verified feature lookups in the caller's atomic transaction.

    Cross-references narrow candidates before reading peptides. Each distinct
    peptide is decompressed at most once from a bounded chunk cache; unmatched
    UniProt isoforms and changed Ensembl sequences are counted and omitted.
    """
    from .data import LOG, fetch_sequence, input_progress, sqlite_activity

    candidates: dict[str, set[str]] = {}
    for stable_id, sequence_id in db.execute(
        "SELECT stable_id, sequence_id FROM sequences WHERE kind='pep'"
    ):
        candidates.setdefault(stable_id, set()).add(sequence_id)
    xrefs: dict[str, set[str]] = {}
    tables = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "ensembl_object_xref" in tables and {"xref_id", "dbprimary_acc"} <= {
        row[1] for row in db.execute("PRAGMA table_info(ensembl_xref)")
    }:
        # Some tiny fixtures do not have a full object_xref table.
        columns = {row[1] for row in db.execute("PRAGMA table_info(ensembl_object_xref)")}
        if {"xref_id", "ensembl_id", "ensembl_object_type"} <= columns:
            external_join = ""
            external_filter = ""
            if (
                "ensembl_external_db" in tables
                and "external_db_id"
                in {row[1] for row in db.execute("PRAGMA table_info(ensembl_xref)")}
                and {"external_db_id", "db_name"}
                <= {row[1] for row in db.execute("PRAGMA table_info(ensembl_external_db)")}
            ):
                external_join = "JOIN ensembl_external_db d ON d.external_db_id=x.external_db_id"
                external_filter = "AND LOWER(d.db_name) LIKE '%uniprot%'"
            with sqlite_activity("Read UniProt protein cross-references"):
                for accession, stable_id in db.execute(f"""
                    SELECT x.dbprimary_acc, t.stable_id FROM ensembl_xref x
                    JOIN ensembl_object_xref o USING (xref_id)
                    JOIN ensembl_translation t ON t.translation_id=o.ensembl_id
                    {external_join}
                    WHERE o.ensembl_object_type='Translation' {external_filter}
                """):
                    if stable_id in candidates:
                        xrefs.setdefault(accession, set()).add(stable_id)
    db.execute("DROP TABLE IF EXISTS ff_uniprot_features")
    db.execute(
        "CREATE TABLE ff_uniprot_features (sequence_id TEXT NOT NULL, feature_id TEXT NOT NULL, "
        "payload TEXT NOT NULL, PRIMARY KEY(sequence_id, feature_id)) WITHOUT ROWID"
    )
    counts = {
        "entries": 0,
        "reviewed_human_entries": 0,
        "matched_proteins": 0,
        "features": 0,
        "sequence_mismatches": 0,
    }
    matched: set[str] = set()
    cache: OrderedDict[tuple[str, str, int], bytes] = OrderedDict()
    with path.open("rb") as source:
        compressed = source.read(2) == b"\x1f\x8b"
    with input_progress(path, "Import reviewed human UniProt", compressed=compressed) as (
        stream,
        bar,
        update,
    ):
        parser = ET.iterparse(stream, events=("start", "end"))
        _, root = next(parser)
        if root.tag != NS + "uniprot":
            raise ValueError("Expected UniProt XML")
        for event, entry in parser:
            if event != "end" or entry.tag != NS + "entry":
                continue
            proteins = entry_proteins(entry)
            reviewed = False
            for protein in proteins:
                reviewed = True
                stable_ids = set(protein["protein_ids"])
                for accession in protein["accessions"]:
                    stable_ids.update(xrefs.get(accession, ()))
                for stable_id in stable_ids:
                    for sequence_id in candidates.get(stable_id, ()):
                        peptide = fetch_sequence(db, "pep", sequence_id, _chunk_cache=cache)
                        if peptide != protein["sequence"]:
                            counts["sequence_mismatches"] += 1
                            continue
                        matched.add(sequence_id)
                        for feature in protein["features"]:
                            payload = json.dumps(feature, separators=(",", ":"))
                            db.execute(
                                "INSERT OR IGNORE INTO ff_uniprot_features VALUES (?, ?, ?)",
                                (sequence_id, feature["feature_id"], payload),
                            )
            counts["entries"] += 1
            counts["reviewed_human_entries"] += int(reviewed)
            update()
            if counts["entries"] % 100 == 0:
                bar.set_postfix(entries=counts["entries"], matched=len(matched), refresh=False)
            root.clear()  # Avoid retaining the full human XML tree in memory.
    if counts["reviewed_human_entries"] == 0:
        raise ValueError("UniProt XML contains no reviewed human entries")
    counts["matched_proteins"] = len(matched)
    counts["features"] = db.execute("SELECT COUNT(*) FROM ff_uniprot_features").fetchone()[0]
    LOG.info("UniProt import: %s", counts)
    return counts
