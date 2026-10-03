"""Compact reference exports and checksum-pinned, resumable installation.

The publication catalog contains exact URLs, not a Zenodo username. Only
preparation uses it; runtime annotation continues to use local SQLite files.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import sqlite3
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Mapping, Sequence
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import data

if TYPE_CHECKING:
    from http.client import HTTPResponse

CATALOG_URL = "https://raw.githubusercontent.com/creisle/fusion_function/main/src/fusion_function/reference_catalog.json"
CATALOG_PATH = Path(__file__).with_name("reference_catalog.json")
RUNTIME_TABLES = ("build_metadata", "ff_transcripts", "ff_interpro", "sequences", "sequence_chunks")
COMPATIBILITY = {
    "format_version": "1",
    "preprocessing_version": "1",
    "cds_mapping_version": "2",
    "feature_annotation_version": "3",
    "transcript_payload_codec": data.TRANSCRIPT_PAYLOAD_CODEC,
    "sequence_chunk_size": str(data.CHUNK_SIZE),
    "sequence_codec": "zlib",
}
DATA_NOTICES = {
    "Ensembl": {
        "terms": "Unrestricted project-generated data; third-party constraints may apply.",
        "url": "https://www.ensembl.org/info/about/legal/disclaimer.html",
    },
    "UniProt Consortium": {
        "terms": "CC BY 4.0. Credit UniProt, link to the license and identify modifications.",
        "url": "https://www.uniprot.org/help/license",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
    },
    "InterPro Consortium": {
        "terms": "Current InterPro downloads: CC0 1.0. Retain historical-source notices.",
        "url": "https://interpro-documentation.readthedocs.io/en/latest/license.html",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
    },
    "PANTHER": {
        "terms": "Classification release 14.1 and 17.0 READMEs carry GPL-2.0-or-later notices. Confirm terms for the imported release and derived classifications before redistribution.",
        "url": "https://data.pantherdb.org/ftp/sequence_classifications/",
    },
    "Ensembl member annotations": {
        "terms": "Member-source terms are not replaced by Ensembl or InterPro terms. PROSITE database terms are CC BY-NC-ND 4.0 with commercial licensing; SMART models require a license. Confirm terms for derived match annotations.",
        "url": "https://prosite.expasy.org/prosite_license.html",
        "smart_url": "https://smart.embl.de/about.cgi",
    },
}


def https_url(value: object) -> str:
    """Require public HTTPS URLs without embedded credentials."""
    if not isinstance(value, str):
        raise ValueError("Reference URL must be an HTTPS URL")
    parsed = urllib.parse.urlparse(value)
    if (
        parsed.scheme != "https"
        or not parsed.netloc
        or parsed.username
        or parsed.password
        or parsed.fragment
    ):
        raise ValueError("Reference URL must be HTTPS without credentials or a fragment")
    return value


def read_catalog(source: str | None = None) -> list[dict[str, Any]]:
    """Refresh the default catalog, falling back to its bundled copy on outages.

    Explicit catalogs and malformed data fail clearly. A failed artifact
    download never silently starts an expensive source build.
    """
    selected = source or os.environ.get("FUSION_FUNCTION_REFERENCE_CATALOG")
    if selected and "://" not in selected:
        raw = Path(selected).expanduser().read_bytes()
    else:
        try:
            with data.open_url(https_url(selected or CATALOG_URL)) as response:
                raw = response.read(4 * data.CHUNK_SIZE + 1)
            if len(raw) > 4 * data.CHUNK_SIZE:
                raise ValueError("Reference catalog exceeds 4 MiB")
        except (OSError, urllib.error.URLError) as exc:
            if selected:
                raise
            data.LOG.warning("Reference catalog unavailable (%s); using bundled catalog", exc)
            raw = CATALOG_PATH.read_bytes()
    catalog = json.loads(raw)
    if (
        not isinstance(catalog, dict)
        or catalog.get("catalog_version") != 1
        or not isinstance(catalog.get("references"), list)
    ):
        raise ValueError(
            "Invalid reference catalog; expected catalog_version=1 and references list"
        )
    entries = []
    for entry in catalog["references"]:
        if not isinstance(entry, dict):
            raise ValueError("Reference catalog entries must be objects")
        # Export manifests can exist before their Zenodo record is published.
        if entry.get("url") is None:
            continue
        https_url(entry["url"])
        for key in ("release", "build_revision", "archive_bytes", "database_bytes"):
            if type(entry.get(key)) is not int or entry[key] <= 0:
                raise ValueError(f"Reference catalog {key} must be a positive integer")
        for key in ("archive_sha256", "database_sha256"):
            if not isinstance(entry.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", entry[key]):
                raise ValueError(f"Invalid reference catalog {key}")
        if entry.get("compression") != "gzip" or not isinstance(entry.get("metadata"), dict):
            raise ValueError("Reference catalog requires gzip compression and metadata")
        entries.append(entry)
    return entries


def select_reference(
    entries: Sequence[Mapping[str, Any]], release: int
) -> Mapping[str, Any] | None:
    """Choose the newest compatible build of the exact requested release."""
    candidates = [
        entry
        for entry in entries
        if entry["release"] == release
        and entry.get("species") == data.SPECIES
        and entry.get("assembly") == data.ASSEMBLY
        and all(entry["metadata"].get(key) == value for key, value in COMPATIBILITY.items())
    ]
    if not candidates:
        return None
    revision = max(entry["build_revision"] for entry in candidates)
    newest = [entry for entry in candidates if entry["build_revision"] == revision]
    if len(newest) != 1:
        raise ValueError(f"Ambiguous reference for release {release}, build {revision}")
    return newest[0]


def open_download(request: urllib.request.Request) -> HTTPResponse:
    """Open an anonymous download, including Range headers for interrupted files."""
    return urllib.request.urlopen(request, timeout=60)  # type: ignore[no-any-return]


def download_archive(entry: Mapping[str, Any], cache: Path) -> Path:
    """Lock the shared archive cache, so installs to different outputs cannot race."""
    directory = cache / "prebuilt" / entry["archive_sha256"]
    directory.mkdir(parents=True, exist_ok=True)
    with data.build_lock(directory / ".download.lock"):
        return _download_archive(entry, directory)


def _download_archive(entry: Mapping[str, Any], directory: Path) -> Path:
    """Resume partial bytes; restart when Range is ignored and verify checksums.

    Interrupted transfers retain partial bytes. Complete but invalid files are
    deleted so the next attempt can start again.
    """
    archive = directory / "ensembl.sqlite.gz"
    partial = archive.with_name(archive.name + ".part")
    size = entry["archive_bytes"]
    data.LOG.info("Prebuilt archive: %s; resumable partial: %s", archive, partial)
    if not archive.exists():
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > size:
            partial.unlink()
            offset = 0
        if offset < size:
            request = urllib.request.Request(
                entry["url"], headers={"Range": f"bytes={offset}-"} if offset else {}
            )
            with open_download(request) as response:
                if response.status == 206:
                    if response.headers.get("Content-Range") != f"bytes {offset}-{size - 1}/{size}":
                        raise ValueError("Unexpected Content-Range for prebuilt reference")
                elif response.status == 200:
                    offset = 0
                else:
                    raise OSError(f"Unexpected reference download status: {response.status}")
                with (
                    partial.open("ab" if offset else "wb") as stream,
                    data.progress(
                        total=size, desc="Download prebuilt reference", unit="B", unit_scale=True
                    ) as bar,
                ):
                    bar.update(offset)
                    received = offset
                    while block := response.read(4 * data.CHUNK_SIZE):
                        if received + len(block) > size:
                            raise ValueError("Reference archive exceeds the catalog size")
                        stream.write(block)
                        received += len(block)
                        bar.update(len(block))
                if received != size:
                    raise OSError(
                        f"Incomplete reference download: {received}/{size} bytes; rerun to resume"
                    )
        if (
            partial.stat().st_size != size
            or data.sha256_file(partial, show_progress=True) != entry["archive_sha256"]
        ):
            partial.unlink()
            raise ValueError("Prebuilt archive checksum/size mismatch; rerun to download again")
        partial.replace(archive)
    elif (
        archive.stat().st_size != size
        or data.sha256_file(archive, show_progress=True) != entry["archive_sha256"]
    ):
        archive.unlink()
        raise ValueError("Cached prebuilt archive checksum/size mismatch; rerun to download again")
    return archive


def validate_database(path: Path, release: int, metadata: Mapping[str, str]) -> None:
    """Check metadata and runtime tables without another expensive integrity scan."""
    from .reference import ReferenceDatabase

    with ReferenceDatabase(path, release=release) as reader:
        for key, value in metadata.items():
            if reader.metadata.get(key) != value:
                raise ValueError(f"Prebuilt database metadata mismatch: {key}")
        tables = {
            row[0] for row in reader.db.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }
        if set(RUNTIME_TABLES) - tables:
            raise ValueError("Prebuilt database is missing runtime tables")
        if (
            reader.db.execute(
                "SELECT 1 FROM ff_transcripts WHERE status='ready' LIMIT 1"
            ).fetchone()
            is None
        ):
            raise ValueError("Prebuilt database contains no usable coding transcripts")
        if (
            reader.db.execute("SELECT 1 FROM sequence_chunks WHERE kind='dna' LIMIT 1").fetchone()
            is None
        ):
            raise ValueError("Prebuilt database contains no genome sequence chunks")


def install_reference(
    *,
    release: int | None,
    cache_dir: Path,
    output: Path | None = None,
    force: bool = False,
    catalog: str | None = None,
) -> Path | None:
    """Install a published reference, or return None when a source build is needed."""
    entries = read_catalog(catalog)
    if not entries:
        data.LOG.info("No published prebuilt references in the catalog; building from Ensembl FTP")
        return None
    if release is None:
        # Latest still means latest Ensembl, not latest uploaded database.
        release, _, _ = data.resolve_core(data.BASE_URL, None)
    entry = select_reference(entries, release)
    if entry is None:
        data.LOG.info("No compatible prebuilt for Ensembl %s; building from Ensembl FTP", release)
        return None
    path = (
        (
            output
            or cache_dir / data.SPECIES / data.ASSEMBLY / f"release-{release}" / "ensembl.sqlite"
        )
        .expanduser()
        .resolve()
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    data.LOG.info(
        "Selected prebuilt: Ensembl %s, build %s, %s",
        release,
        entry["build_revision"],
        entry["url"],
    )
    data.LOG.info(
        "Final reference: %s; expected installed size: %s bytes", path, entry["database_bytes"]
    )
    staging = path.with_name(path.name + ".prebuilt.part")
    data.LOG.info("Expanded partial database: %s", staging)
    steps = data.BuildSteps(4)
    with data.build_lock(path.with_name(path.name + ".prepare.lock")):
        if path.exists() and not force:
            validate_database(path, release, COMPATIBILITY)
            data.LOG.info("Using existing reference: %s; --force replaces it", path)
            return path
        try:
            with steps.step("Download and verify prebuilt archive"):
                data.LOG.info("A large reference download can take over 10 minutes")
                archive = download_archive(entry, cache_dir)
            with steps.step("Expand reference database"):
                data.LOG.info("Expanding a large reference can take over 10 minutes")
                digest = hashlib.sha256()
                size = 0
                with (
                    gzip.open(archive, "rb") as source,
                    staging.open("wb") as target,
                    data.progress(
                        total=entry["database_bytes"],
                        desc="Expand prebuilt reference",
                        unit="B",
                        unit_scale=True,
                    ) as bar,
                ):
                    while block := source.read(4 * data.CHUNK_SIZE):
                        size += len(block)
                        if size > entry["database_bytes"]:
                            raise ValueError("Expanded reference exceeds the catalog size")
                        digest.update(block)
                        target.write(block)
                        bar.update(len(block))
                if (
                    size != entry["database_bytes"]
                    or digest.hexdigest() != entry["database_sha256"]
                ):
                    raise ValueError("Expanded reference checksum/size mismatch")
            with steps.step("Validate reference format and release"):
                validate_database(staging, release, entry["metadata"])
            with steps.step("Install completed reference"):
                os.replace(staging, path)
        finally:
            staging.unlink(missing_ok=True)
    steps.complete()
    return path


def export_reference(
    database: Path, directory: Path, *, build_revision: int = 1, zenodo_record: int | None = None
) -> Path:
    """Export a read-only snapshot of runtime data; preserve the maintainer's DB.

    Bulk INSERT SELECT copies compressed payloads without Python decoding. The
    attached source is read-only and one transaction holds a consistent snapshot.
    """
    if build_revision <= 0 or (zenodo_record is not None and zenodo_record <= 0):
        raise ValueError("Build revision and Zenodo record must be positive")
    database = database.expanduser().resolve()
    directory = directory.expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    with data.ReferenceReader(database) as reader:
        metadata = dict(reader.metadata)
    release = int(metadata["release"])
    validate_database(database, release, COMPATIBILITY)
    stem = f"{data.SPECIES}.{data.ASSEMBLY}.ensembl-{release}.ff-v1.build-{build_revision}"
    archive, manifest = directory / (stem + ".sqlite.gz"), directory / (stem + ".json")
    with data.build_lock(directory / (stem + ".export.lock")):
        if archive.exists() or manifest.exists():
            raise FileExistsError("Export already exists; use a new build revision or directory")
        with tempfile.TemporaryDirectory(prefix=stem + ".", dir=directory) as workspace:
            compact = Path(workspace) / "ensembl.sqlite"
            data.LOG.info("Export destination: %s; temporary database: %s", archive, compact)
            steps = data.BuildSteps(4)
            with closing(sqlite3.connect(compact, uri=True)) as target:
                target.execute("ATTACH DATABASE ? AS original", (database.as_uri() + "?mode=ro",))
                target.execute("PRAGMA secure_delete=OFF")
                target.execute("PRAGMA user_version=1")
                target.execute("BEGIN")
                with steps.step("Copy runtime tables into a compact database"):
                    data.LOG.info("Copying a large reference can take over 10 minutes")
                    for table in (*RUNTIME_TABLES, "source_files"):
                        schema = target.execute(
                            "SELECT sql FROM original.sqlite_master WHERE type='table' AND name=?",
                            (table,),
                        ).fetchone()
                        if schema is None:
                            if table == "source_files":
                                continue
                            raise ValueError(f"Missing runtime table: {table}")
                        target.execute(schema[0])
                        data.sqlite_phase(
                            target,
                            f'INSERT INTO "{table}" SELECT * FROM original."{table}"',
                            f"Export {table}",
                        )
                        for (sql,) in target.execute(
                            "SELECT sql FROM original.sqlite_master WHERE type='index' AND tbl_name=? AND sql IS NOT NULL",
                            (table,),
                        ).fetchall():
                            data.sqlite_phase(target, sql, f"Index exported {table}")
                    # Keep original provenance, plus notices describing our transformation.
                    metadata = dict(target.execute("SELECT * FROM build_metadata"))
                    target.executemany(
                        "INSERT OR REPLACE INTO build_metadata VALUES (?,?)",
                        [
                            ("reference_kind", "runtime"),
                            ("reference_build_revision", str(build_revision)),
                            ("exported_utc", datetime.now(timezone.utc).isoformat()),
                            ("data_source_notices", json.dumps(DATA_NOTICES, sort_keys=True)),
                            (
                                "data_modifications",
                                "Runtime subset; derived transcript, feature and splice mappings; lossless compression.",
                            ),
                        ],
                    )
                    target.commit()
                with steps.step("Validate exported SQLite integrity"):
                    if data.sqlite_phase(target, "PRAGMA integrity_check", "Check database") != [
                        ("ok",)
                    ]:
                        raise ValueError("Exported SQLite integrity check failed")
            validate_database(compact, release, COMPATIBILITY)
            with steps.step("Compress exported reference"):
                data.LOG.info("Compressing a large reference can take over 10 minutes")
                digest = hashlib.sha256()
                temporary_archive = Path(workspace) / archive.name
                with (
                    compact.open("rb") as source,
                    temporary_archive.open("wb") as output,
                    gzip.GzipFile(
                        filename="", fileobj=output, mode="wb", compresslevel=1, mtime=0
                    ) as compressed,
                    data.progress(
                        total=compact.stat().st_size,
                        desc="Compress reference",
                        unit="B",
                        unit_scale=True,
                    ) as bar,
                ):
                    while block := source.read(4 * data.CHUNK_SIZE):
                        digest.update(block)
                        compressed.write(block)
                        bar.update(len(block))
            with steps.step("Write manifest and publish export"):
                entry = {
                    "species": data.SPECIES,
                    "assembly": data.ASSEMBLY,
                    "release": release,
                    "build_revision": build_revision,
                    "url": f"https://zenodo.org/records/{zenodo_record}/files/{archive.name}?download=1"
                    if zenodo_record
                    else None,
                    "filename": archive.name,
                    "compression": "gzip",
                    "archive_bytes": temporary_archive.stat().st_size,
                    "archive_sha256": data.sha256_file(temporary_archive, show_progress=True),
                    "database_bytes": compact.stat().st_size,
                    "database_sha256": digest.hexdigest(),
                    "metadata": COMPATIBILITY,
                    "data_source_notices": DATA_NOTICES,
                    "source_metadata": metadata,
                }
                temporary_manifest = Path(workspace) / manifest.name
                temporary_manifest.write_text(
                    json.dumps({"catalog_version": 1, "references": [entry]}, indent=2) + "\n"
                )
                temporary_archive.replace(archive)
                temporary_manifest.replace(manifest)
            steps.complete()
    if zenodo_record is None:
        data.LOG.info(
            "Manifest URL is unset; add the published Zenodo file URL before registering it"
        )
    data.LOG.info("Export archive: %s; manifest: %s", archive, manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    """Maintainer CLI to export a tested reference for publication."""
    from .reference import ReferenceDatabase

    parser = argparse.ArgumentParser(prog="fusion-function export-reference")
    parser.add_argument(
        "database",
        type=Path,
        nargs="?",
        help="Existing database; default: FUSION_FUNCTION_DB or newest installed local release",
    )
    parser.add_argument(
        "--release", type=int, help="Installed Ensembl release; default: newest installed release"
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Directory for compressed DB and manifest"
    )
    parser.add_argument("--build-revision", type=int, default=1)
    parser.add_argument(
        "--zenodo-record", type=int, help="Reserved/published Zenodo record ID; not a username"
    )
    args = parser.parse_args(argv)
    if args.release is not None and args.release <= 0:
        parser.error("--release must be positive")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        # Reuse analysis selection, including cache overrides and release checks.
        # Selection reads local files only; export never builds or downloads a DB.
        with ReferenceDatabase(args.database, release=args.release) as reference:
            database = reference.path
        data.LOG.info("Export source database: %s", database)
        manifest = export_reference(
            database,
            args.output,
            build_revision=args.build_revision,
            zenodo_record=args.zenodo_record,
        )
    except (OSError, ValueError, sqlite3.Error) as exc:
        data.LOG.error("Export failed: %s", exc)
        return 1
    except KeyboardInterrupt:
        data.LOG.error("Export interrupted; source database unchanged")
        return 130
    print(manifest)
    return 0
