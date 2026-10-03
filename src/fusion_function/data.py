#!/usr/bin/env python3
"""Build and read the human GRCh38 reference used by fusion-function.

Commands (Python 3.11+, tqdm is the only external dependency):
    fusion-function prepare-data                      # latest Ensembl release
    fusion-function prepare-data --release 116        # pinned release
    fusion-function prepare-data --preprocess-only DB # regenerate derived data

Default output: <cache>/homo_sapiens/GRCh38/release-<release>/ensembl.sqlite.
Set FUSION_FUNCTION_CACHEDIR or --cache-dir to change the cache root.

Storage: ensembl_* tables preserve the FTP dump columns. DNA and peptides use
independently compressed 1 MiB chunks. ff_transcripts stores prepared exon/CDS
coordinates, splice windows and protein features; ff_interpro stores names and
types. Transcript models use compressed JSON; error messages remain plain JSON.
ff_uniprot_features stores sequence-verified reviewed human sites,
domains and motifs. Coordinates are 1-based and inclusive; CDS blocks follow transcript
orientation, including on the negative strand. Pre-mRNA is reconstructed when
requested, rather than stored for every transcript.

Preparation downloads FTP files over HTTPS and verifies cached files by SHA-256.
InterPro metadata is fetched by default, with archived entry lists supplementing
missing accessions. --interpro-entries FILE supplies local metadata instead.
Archive checksums and source releases are recorded in build_metadata.
Human PANTHER subfamily classifications are fetched for the version recorded
by Ensembl, or supplied with --panther-classifications FILE. Protein cross-
references enrich existing family hits; no domain coordinates are invented.
Reviewed human UniProt XML is fetched during preparation only, or supplied with
--uniprot-features FILE. Features require exact protein sequence identity and
carry UniProt accessions, isoforms and evidence codes.

New databases resume from <output>.building and are published only after
integrity checking; --force permits replacement. Repeating the same build skips
completed imports and reuses committed FASTA records from unchanged inputs.
Reprocessing updates derived tables in one transaction and rolls
back on failure. Numbered steps and terminal progress bars show build progress.
Runtime readers are local and read-only; no REST API, MySQL server or pysam is
required. See README.md and docs/reference.md for configuration details.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import re
import sqlite3
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from bisect import bisect_left
from collections import OrderedDict
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import AbstractContextManager, ExitStack, closing, contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from html.parser import HTMLParser
from itertools import groupby
from pathlib import Path
from types import TracebackType
from typing import TYPE_CHECKING, Any, Literal, Self, TextIO, TypedDict, cast, overload
from xml.etree.ElementTree import ParseError

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

if TYPE_CHECKING:
    from http.client import HTTPResponse

    from .ensembl import (
        CDSBlock,
        EnsemblError,
        GenomicSegment,
        ProteinFeatureAnnotation,
        ProteinFeatureResponse,
        TranscriptProteinFeatureResult,
    )


# Raw SQL rows have different columns per source table; derived structures use
# the same typed records as the annotation API.
DatabaseRow = dict[str, Any]
AnalysisMetadata = dict[int, tuple[str | None, str | None]]
InterProMetadata = dict[str, tuple[str | None, str | None]]
ChunkCache = OrderedDict[tuple[str, str, int], bytes]
STRUCTURE_SOURCES = frozenset({"sifts", "alphafold"})
TRANSCRIPT_PAYLOAD_CODEC = "zlib-json-v1"


class SourceRecord(TypedDict):
    url: str
    sha256: str
    bytes: int


class TranscriptErrorSummary(TypedDict):
    count: int
    example_transcripts: list[str]


class CheckpointRecord(TypedDict):
    source: dict[str, Any]
    fingerprint: str
    count: int
    complete: bool


BASE_URL = "https://ftp.ensembl.org/pub/"
INTERPRO_ENTRIES_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
INTERPRO_RELEASES_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/"
PANTHER_CLASSIFICATIONS_URL = "https://data.pantherdb.org/ftp/sequence_classifications/"
# These inputs are fixed: omitting a table or sequence type can break preprocessing.
TABLES = (
    "meta",
    "coord_system",
    "seq_region",
    "gene",
    "transcript",
    "exon",
    "exon_transcript",
    "translation",
    "analysis",
    "protein_feature",
    "interpro",
    "xref",
    "external_db",
    "object_xref",
)
PREPROCESS_TABLES = (
    "transcript",
    "translation",
    "exon",
    "exon_transcript",
    "seq_region",
    "coord_system",
    "protein_feature",
    "analysis",
    "interpro",
    "xref",
)
SEQUENCE_TYPES = ("dna", "pep")
SPECIES = "homo_sapiens"
ASSEMBLY = "GRCh38"
CHUNK_SIZE = 1024 * 1024
LOG = logging.getLogger("ensembl_to_sqlite")
STEP_PREFIX = ContextVar("fusion_function_prepare_step", default="")


# Progress reporting and cache configuration.


class BuildSteps:
    """Label the fixed build stages; their durations differ, so no total ETA."""

    def __init__(self, total: int) -> None:
        """Track the expected stage count and total elapsed build time."""
        self.total = total
        self.current = 0
        self.started = time.monotonic()

    @contextmanager
    def step(self, description: str) -> Iterator[None]:
        """Label a stage and its nested bars; restore the label even on failure."""
        self.current += 1
        if self.current > self.total:
            raise ValueError("Preparation step count exceeds its declared total")
        LOG.info("[Step %s/%s] %s", self.current, self.total, description)
        token = STEP_PREFIX.set(f"[{self.current}/{self.total}] ")
        try:
            yield
        finally:
            STEP_PREFIX.reset(token)

    def complete(self) -> None:
        """Check that every declared stage ran, then log the total duration."""
        if self.current != self.total:
            raise ValueError("Preparation completed fewer steps than expected")
        LOG.info(
            "Completed all %s steps in %s",
            self.total,
            tqdm.format_interval(time.monotonic() - self.started),
        )


def progress(
    iterable: Iterable[Any] | None = None,
    *,
    total: int | None = None,
    desc: str = "",
    unit: str = "it",
    unit_scale: bool = False,
    unit_divisor: int = 1000,
    disable: bool | None = None,
    file: TextIO | None = None,
) -> tqdm:
    """Use compact, rate-limited bars; suppress terminal redraws in log files."""
    long_process_notice(desc)
    return tqdm(
        iterable,
        total=total,
        desc=STEP_PREFIX.get() + desc,
        unit=unit,
        unit_scale=unit_scale,
        unit_divisor=unit_divisor,
        disable=disable,
        file=file,
        dynamic_ncols=True,
        mininterval=0.5,
    )


def file_label(path: Path) -> str:
    """Short labels keep progress output on one terminal line."""
    if path.name.endswith(".sql.gz"):
        return "core schema"
    if path.name == "stream":
        return "UniProt features"
    if path.name == "entry.list":
        return "InterPro"
    for kind in SEQUENCE_TYPES:
        if f".{kind}." in path.name:
            return f"{kind} FASTA"
    return path.name.removesuffix(".txt.gz")[:24]


@overload
def input_progress(
    path: Path, description: str, *, compressed: bool = True, text: Literal[True]
) -> AbstractContextManager[tuple[io.TextIOWrapper, tqdm, Callable[[], None]]]: ...


@overload
def input_progress(
    path: Path, description: str, *, compressed: bool = True, text: Literal[False] = False
) -> AbstractContextManager[tuple[gzip.GzipFile | io.BufferedReader, tqdm, Callable[[], None]]]: ...


@overload
def input_progress(
    path: Path, description: str, *, compressed: bool = True, text: bool
) -> AbstractContextManager[
    tuple[io.TextIOWrapper | gzip.GzipFile | io.BufferedReader, tqdm, Callable[[], None]]
]: ...


@contextmanager
def input_progress(
    path: Path, description: str, *, compressed: bool = True, text: bool = False
) -> Iterator[tuple[Any, tqdm, Callable[[], None]]]:
    """Track source bytes during parsing, with no pre-count or gzip size guess."""
    with ExitStack() as stack:
        raw = stack.enter_context(path.open("rb"))
        source = stack.enter_context(gzip.GzipFile(fileobj=raw)) if compressed else raw
        stream = (
            stack.enter_context(io.TextIOWrapper(source, encoding="utf-8", newline=""))
            if text
            else source
        )
        bar = stack.enter_context(
            progress(
                total=path.stat().st_size,
                desc=description,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )
        )
        position = 0

        def update() -> None:
            """Advance by source bytes read, including gzip buffering."""
            nonlocal position
            consumed = raw.tell()
            bar.update(consumed - position)
            position = consumed

        yield stream, bar, update
        update()


def long_process_notice(description: str) -> None:
    """Print one advance notice for work that can exceed ten minutes.

    These are possibilities for large human references, not completion estimates.
    Emit notices only when work starts, so cached or checkpointed stages stay quiet.
    """
    if description == "Drop previous ff_transcripts":
        note = "Removing old transcript models can take an hour or longer on some systems."
    elif description == "Preprocess transcripts":
        note = "Preparing all human transcripts can take an hour or longer."
    elif description == "Check database":
        note = "The full database integrity check can take over 10 minutes."
    elif (
        description.startswith(("Download ", "Import ", "Index "))
        and description not in {"Import InterPro", "Import PANTHER", "Index analysis", "Index meta"}
    ) or description in {
        "Analyze database",
        "Read UniProt protein cross-references",
        "Read genome and protein sequence lengths",
        "Prepare ordered exon stream",
        "Prepare ordered protein feature stream",
        "Prepare transcript stream",
        "Commit preprocessed reference",
        "Roll back preprocessing transaction",
    }:
        note = (
            "This operation can take over 10 minutes, depending on reference size and system load."
        )
    else:
        return
    LOG.info("%s%s: %s", STEP_PREFIX.get(), description, note)


@contextmanager
def sqlite_activity(description: str) -> Iterator[None]:
    """Log phase boundaries without background output or periodic redraws."""
    label = STEP_PREFIX.get() + description
    LOG.info("%s", label)
    long_process_notice(description)
    outcome = "failed"
    try:
        yield
        outcome = "done"
    except KeyboardInterrupt:
        outcome = "interrupted"
        raise
    finally:
        LOG.info("%s: %s", label, outcome)


def sqlite_phase(db: sqlite3.Connection, sql: str, description: str) -> list[tuple[Any, ...]]:
    """Run one SQLite statement with start/completion logs and signal handling."""
    with sqlite_activity(description):

        def checkpoint() -> int:
            """Keep Python signal handling responsive during long SQLite statements."""
            return 0

        db.set_progress_handler(checkpoint, 100000)
        try:
            return db.execute(sql).fetchall()
        finally:
            db.set_progress_handler(None, 0)


def default_cache_dir() -> Path:
    """Use the package override, then the platform's conventional cache root."""
    if value := os.environ.get("FUSION_FUNCTION_CACHEDIR"):
        return Path(value).expanduser()
    if sys.platform == "win32":
        root = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local"))
    elif sys.platform == "darwin":
        root = Path.home() / "Library/Caches"
    else:
        root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return root / "fusion_function"


# Source discovery, verified downloads and MySQL dump imports.


def open_url(url: str) -> HTTPResponse:
    """Open a source URL with an identifying user agent and a bounded timeout."""
    request = urllib.request.Request(url, headers={"User-Agent": "ensembl-to-sqlite/1"})
    return urllib.request.urlopen(request, timeout=60)


class Links(HTMLParser):
    """Collect filenames from the HTML directory listings used by FTP mirrors."""

    def __init__(self) -> None:
        """Start with no names; repeated links are collapsed into a set."""
        super().__init__()
        self.names: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        """Keep link basenames, excluding navigation to the current/parent path."""
        if tag == "a":
            href = dict(attrs).get("href", "")
            path = urllib.parse.unquote(urllib.parse.urlparse(href).path)
            name = path.rstrip("/").rsplit("/", 1)[-1]
            if name not in {"", ".", ".."}:
                self.names.add(name)


def listing(url: str) -> set[str]:
    """Return directory entry names without depending on HTML table formatting."""
    with open_url(url) as response:
        text = response.read().decode("utf-8")
    parser = Links()
    parser.feed(text)
    return parser.names


def resolve_core(base: str, release: int | None) -> tuple[int, str, str]:
    """Find a published release and its human GRCh38 core directory.

    Latest means the largest numbered release, not an unversioned FTP alias.
    Choosing the exact assembly suffix avoids mixing GRCh37 and GRCh38 data.
    """
    if release is None:
        releases = [
            int(m[1]) for name in listing(base) if (m := re.fullmatch(r"release-(\d+)", name))
        ]
        if not releases:
            raise ValueError("No numbered Ensembl releases found; supply --release")
        release = max(releases)
    root = f"{base}release-{release}/"
    core_name = f"{SPECIES}_core_{release}_38"
    if core_name not in listing(root + "mysql/"):
        raise ValueError(f"Human GRCh38 core database {core_name} is absent from release {release}")
    return release, root, core_name


def sha256_file(path: Path, *, show_progress: bool = False) -> str:
    """Hash a file in bounded reads, optionally displaying verification progress."""
    digest = hashlib.sha256()
    with (
        path.open("rb") as stream,
        progress(
            total=path.stat().st_size,
            desc="Verify " + file_label(path),
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            disable=None if show_progress else True,
        ) as bar,
    ):
        while block := stream.read(4 * CHUNK_SIZE):
            digest.update(block)
            bar.update(len(block))
    return digest.hexdigest()


def download(url: str, directory: Path) -> tuple[Path, str]:
    """Reuse verified cached bytes or atomically publish a complete download.

    The caller supplies a directory that separates releases/source locations.
    Recoverable failures are retried; partial files are removed on every exit.
    """
    directory = directory.expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    name = urllib.parse.urlparse(url).path.rsplit("/", 1)[-1]
    path = directory / name
    sidecar = path.with_name(path.name + ".sha256")
    if path.is_file() and sidecar.is_file():
        cached_digest = sha256_file(path, show_progress=True)
        if cached_digest == sidecar.read_text().strip():
            LOG.info("Using cached %s", path)
            return path, cached_digest
        LOG.warning("Cached checksum mismatch: downloading %s again", name)
    LOG.info("Downloading %s", url)
    LOG.info("  Cached file: %s; checksum: %s", path, sidecar)
    for attempt in range(3):
        temporary = None
        try:
            digest = hashlib.sha256()
            received = 0
            with open_url(url) as response:
                expected = response.headers.get("Content-Length")
                with tempfile.NamedTemporaryFile(
                    dir=directory, prefix=name + ".", suffix=".part", delete=False
                ) as out:
                    temporary = Path(out.name)
                    LOG.info("  Partial download: %s", temporary)
                    with progress(
                        total=int(expected) if expected is not None else None,
                        desc="Download " + file_label(path),
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                    ) as bar:
                        while block := response.read(4 * CHUNK_SIZE):
                            out.write(block)
                            digest.update(block)
                            received += len(block)
                            bar.update(len(block))
            if expected is not None and received != int(expected):
                raise OSError(f"Truncated download: expected {expected}, received {received}")
            temporary.replace(path)
            checksum = digest.hexdigest()
            sidecar.write_text(checksum + "\n")
            return path, checksum
        except (OSError, urllib.error.URLError) as exc:
            if isinstance(exc, urllib.error.HTTPError) and exc.code in {400, 403, 404}:
                raise
            if attempt == 2:
                raise
            LOG.warning("Download failed (%s); retrying", exc)
            time.sleep(2**attempt)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise AssertionError("Unreachable")


def sql_name(value: str) -> str:
    """Quote a SQLite identifier; bound parameters can only represent values."""
    return '"' + value.replace('"', '""') + '"'


def parse_schema(path: Path) -> dict[str, list[tuple[str, str]]]:
    """Read column order/types from MySQL DDL without executing source SQL."""
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        text = stream.read()
    tables = {}
    for match in re.finditer(
        r"CREATE TABLE(?: IF NOT EXISTS)? `([^`]+)`\s*\((.*?)\n\)", text, re.S
    ):
        columns = []
        for column in re.finditer(r"^\s*`([^`]+)`\s+([A-Za-z]+)", match[2], re.M):
            mysql_type = column[2].lower()
            # Ignore MySQL constraints/options; imports need only SQLite affinities.
            if mysql_type in {"tinyint", "smallint", "mediumint", "int", "bigint"}:
                affinity = "INTEGER"
            elif mysql_type in {"float", "double", "decimal"}:
                affinity = "REAL"
            else:
                affinity = "TEXT"
            columns.append((column[1], affinity))
        if columns:
            tables[match[1]] = columns
    if not tables:
        raise ValueError("Cannot parse any CREATE TABLE definitions from the source schema")
    return tables


MYSQL_ESCAPES = {"0": "\0", "b": "\b", "n": "\n", "r": "\r", "t": "\t", "Z": "\x1a", "\\": "\\"}


def mysql_value(value: str) -> str | None:
    r"""Decode MySQL dump escapes, distinguishing SQL NULL (\N) from text."""
    if value == r"\N":
        return None
    return re.sub(r"\\(.)", lambda m: MYSQL_ESCAPES.get(m[1], m[1]), value)


def import_table(
    db: sqlite3.Connection, table: str, columns: Sequence[tuple[str, str]], path: Path
) -> int:
    """Import one gzipped tab-separated dump, add lookup indexes and commit."""
    destination = "ensembl_" + table
    definitions = ", ".join(f"{sql_name(name)} {kind}" for name, kind in columns)
    db.execute(f"CREATE TABLE {sql_name(destination)} ({definitions})")
    statement = f"INSERT INTO {sql_name(destination)} VALUES ({','.join('?' for _ in columns)})"
    count = 0
    batch = []
    with input_progress(path, "Import " + table, text=True) as (stream, bar, update):
        for number, line in enumerate(stream, 1):
            values = line.rstrip("\n").removesuffix("\r").split("\t")
            if len(values) != len(columns):
                raise ValueError(
                    f"{path.name}:{number}: {len(values)} fields; expected {len(columns)}"
                )
            batch.append(tuple(mysql_value(value) for value in values))
            if len(batch) == 10000:
                db.executemany(statement, batch)
                count += len(batch)
                batch.clear()
                bar.set_postfix(rows=f"{count:,}", refresh=False)
                update()
        if batch:
            db.executemany(statement, batch)
            count += len(batch)
        bar.set_postfix(rows=f"{count:,}", refresh=False)
    names = {name for name, _ in columns}
    # Primary lookup IDs plus joins commonly needed for transcript/domain data.
    index_columns = list(
        dict.fromkeys(
            name
            for name in (
                table + "_id",
                "stable_id",
                "transcript_id",
                "translation_id",
                "exon_id",
                "gene_id",
                "seq_region_id",
                "analysis_id",
                "xref_id",
                "ensembl_id",
                "dbprimary_acc",
                "interpro_ac",
                "id",
                "meta_key",
            )
            if name in names
        )
    )
    for name in progress(index_columns, desc="Index " + table, unit="index"):
        db.execute(
            f"CREATE INDEX IF NOT EXISTS {sql_name(destination + '_' + name)} "
            f"ON {sql_name(destination)} ({sql_name(name)})"
        )
    db.commit()
    LOG.info("Imported %s: %s rows", table, f"{count:,}")
    return count


# Chunked sequence storage and inclusive interval reads.


def create_sequence_tables(db: sqlite3.Connection) -> None:
    """Create sequence metadata and separately addressable compressed chunks."""
    db.executescript("""
        CREATE TABLE sequences (
            kind TEXT NOT NULL, sequence_id TEXT NOT NULL, stable_id TEXT NOT NULL,
            header TEXT NOT NULL, length INTEGER NOT NULL, sha256 TEXT NOT NULL,
            PRIMARY KEY (kind, sequence_id)
        ) WITHOUT ROWID;
        CREATE INDEX sequence_stable_id ON sequences (kind, stable_id);
        CREATE TABLE sequence_chunks (
            kind TEXT NOT NULL, sequence_id TEXT NOT NULL,
            start INTEGER NOT NULL, data BLOB NOT NULL,
            PRIMARY KEY (kind, sequence_id, start)
        ) WITHOUT ROWID;
    """)


def validate_dna_regions(db: sqlite3.Connection) -> None:
    """Names are unique within coordinate systems, not across genome assemblies."""
    db.execute("CREATE INDEX IF NOT EXISTS ensembl_seq_region_name ON ensembl_seq_region (name)")
    unmatched = db.execute(
        "SELECT s.sequence_id, s.length FROM sequences s WHERE s.kind='dna' AND NOT EXISTS ("
        "SELECT 1 FROM ensembl_seq_region r JOIN ensembl_coord_system c USING (coord_system_id) "
        "WHERE r.name=s.sequence_id AND c.version=? AND r.length=s.length) LIMIT 1",
        (ASSEMBLY,),
    ).fetchone()
    if unmatched:
        sequence_id, length = unmatched
        candidates = db.execute(
            "SELECT r.length FROM ensembl_seq_region r "
            "JOIN ensembl_coord_system c USING (coord_system_id) "
            "WHERE r.name=? AND c.version=?",
            (sequence_id, ASSEMBLY),
        ).fetchall()
        lengths = sorted({row[0] for row in candidates})
        raise ValueError(
            f"FASTA/core sequence-region mismatch for {sequence_id}: FASTA length {length}; "
            f"{ASSEMBLY} core lengths {lengths or 'no matching region'}"
        )
    LOG.info("Validated DNA sequence-region lengths against %s core records", ASSEMBLY)


def import_fasta(db: sqlite3.Connection, kind: str, path: Path, *, resume: bool = False) -> int:
    """Stream FASTA, optionally reusing committed records from the same input.

    The builder checks the source checksum before enabling resume. Record rows
    exist only after their chunks are complete; discard any orphaned chunks.
    """
    completed = set()
    if resume:
        completed = {
            row[0] for row in db.execute("SELECT sequence_id FROM sequences WHERE kind=?", (kind,))
        }
        db.execute(
            "DELETE FROM sequence_chunks WHERE kind=? AND NOT EXISTS ("
            "SELECT 1 FROM sequences s WHERE s.kind=sequence_chunks.kind "
            "AND s.sequence_id=sequence_chunks.sequence_id)",
            (kind,),
        )
        db.commit()
        if completed:
            LOG.info(
                "Resuming %s FASTA: reusing %s complete sequences; rereading gzip to the next record",
                kind,
                f"{len(completed):,}",
            )
    sequence_id = header = None
    skip_record = False
    buffer = bytearray()
    offset = count = 0
    digest = hashlib.sha256()

    def flush(size: int) -> None:
        """Store the next chunk; its database start is 1-based, not a byte offset."""
        nonlocal offset
        chunk = bytes(buffer[:size])
        del buffer[:size]
        db.execute(
            "INSERT INTO sequence_chunks VALUES (?, ?, ?, ?)",
            (kind, sequence_id, offset + 1, zlib.compress(chunk, level=1)),
        )
        offset += len(chunk)
        if size == CHUNK_SIZE:
            update()

    def finish() -> None:
        """Flush one FASTA record and record its full length, header and checksum."""
        nonlocal count
        if sequence_id is None:
            return
        if skip_record:
            count += 1
            return
        if buffer:
            flush(len(buffer))
        if offset == 0:
            raise ValueError(f"Empty FASTA sequence: {sequence_id}")
        # Assembly sequence names may contain dots; only peptide version suffixes
        # are removed to support lookups by either stable or versioned protein ID.
        stable_id = sequence_id if kind == "dna" else re.sub(r"\.\d+$", "", sequence_id)
        db.execute(
            "INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?)",
            (kind, sequence_id, stable_id, header, offset, digest.hexdigest()),
        )
        count += 1
        if count % 1000 == 0 or kind == "dna":
            db.commit()
            bar.set_postfix(records=f"{count:,}", refresh=False)
            update()

    with input_progress(path, "Import " + kind + " FASTA") as (stream, bar, update):
        for line in stream:
            if line.startswith(b">"):
                finish()
                header = line[1:].strip().decode("utf-8")
                sequence_id = header.split()[0]
                skip_record = sequence_id in completed
                offset = 0
                digest = hashlib.sha256()
            else:
                if skip_record:
                    continue
                sequence = line.strip().upper()
                if not sequence:
                    continue
                if sequence_id is None:
                    raise ValueError("FASTA sequence precedes its header")
                alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ*" if kind == "pep" else b"ACGTRYSWKMBDHVN"
                if sequence.translate(None, alphabet):
                    raise ValueError(f"Invalid {kind} sequence characters in {sequence_id}")
                digest.update(sequence)
                buffer.extend(sequence)
                while len(buffer) >= CHUNK_SIZE:
                    flush(CHUNK_SIZE)
        finish()
        bar.set_postfix(records=f"{count:,}", refresh=False)
    if not count:
        raise ValueError(f"No FASTA records in {path}")
    db.commit()
    LOG.info("Imported %s: %s sequences", kind, f"{count:,}")
    return count


def fetch_sequence(
    db: sqlite3.Connection,
    kind: str,
    sequence_id: str,
    start: int = 1,
    end: int | None = None,
    strand: int = 1,
    *,
    _chunk_cache: ChunkCache | None = None,
    _cache_limit: int = 64,
) -> str:
    """Read a 1-based inclusive interval, optionally using a bounded chunk cache.

    Unversioned protein IDs are accepted only when they match one sequence.
    Negative-strand nucleotide reads return the reverse complement.
    """
    row = db.execute(
        "SELECT sequence_id, length FROM sequences WHERE kind=? AND sequence_id=?",
        (kind, sequence_id),
    ).fetchone()
    if row is None:
        # Never silently choose among several versions of an unversioned ID.
        rows = db.execute(
            "SELECT sequence_id, length FROM sequences WHERE kind=? AND stable_id=?",
            (kind, sequence_id),
        ).fetchall()
        if len(rows) != 1:
            raise KeyError(f"Missing or ambiguous sequence: {kind}/{sequence_id}")
        row = rows[0]
    sequence_id, length = row
    end = length if end is None else end
    if not 1 <= start <= end <= length:
        raise ValueError(f"Invalid interval {start}-{end}; sequence length is {length}")
    if strand not in {-1, 1} or (kind == "pep" and strand != 1):
        raise ValueError("Strand must be +1/-1 for nucleotides or +1 for peptides")
    # Convert the requested position to the 1-based start of its containing chunk.
    first_chunk = ((start - 1) // CHUNK_SIZE) * CHUNK_SIZE + 1
    pieces = []
    for position in range(first_chunk, end + 1, CHUNK_SIZE):
        key = kind, sequence_id, position
        if _chunk_cache is not None and key in _chunk_cache:
            sequence = _chunk_cache[key]
            _chunk_cache.move_to_end(key)
        else:
            record = db.execute(
                "SELECT data FROM sequence_chunks WHERE kind=? AND sequence_id=? AND start=?", key
            ).fetchone()
            if record is None:
                raise ValueError(f"Missing sequence chunk: {key}")
            sequence = zlib.decompress(record[0])
            if _chunk_cache is not None and _cache_limit:
                _chunk_cache[key] = sequence
                while len(_chunk_cache) > _cache_limit:
                    _chunk_cache.popitem(last=False)
        # Python slices are 0-based/exclusive; database intervals are inclusive.
        pieces.append(sequence[max(0, start - position) : end - position + 1])
    result = b"".join(pieces).decode("ascii")
    if len(result) != end - start + 1:
        raise ValueError("Missing or damaged sequence chunks")
    if strand == -1:
        result = result.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
    return result


# Streaming transcript models and coordinate mapping.


def dict_rows(
    db: sqlite3.Connection, query: str, *, description: str | None = None
) -> Iterator[DatabaseRow]:
    """Yield named rows from a query without changing the connection row factory."""
    # Sorted queries can do substantial work before returning their first row.
    if description is not None:
        with sqlite_activity(description):
            cursor = db.execute(query)
    else:
        cursor = db.execute(query)
    names = [column[0] for column in cursor.description]
    for row in cursor:
        yield dict(zip(names, row))


class GroupStream:
    """Consume a query grouped by internal transcript ID without loading it all.

    Input rows and calls to take() must both follow ascending transcript ID.
    This lets exon/feature queries advance alongside the main transcript query.
    """

    def __init__(self, rows: Iterable[DatabaseRow]) -> None:
        """Prime the first group from an already ordered row iterator."""
        self.groups = iter(groupby(rows, key=lambda row: row["transcript_id"]))
        self.current = next(self.groups, None)

    def take(self, transcript_id: int) -> list[DatabaseRow]:
        """Return this transcript's rows, or an empty list when it has no group."""
        while self.current is not None and self.current[0] < transcript_id:
            self.current = next(self.groups, None)
        if self.current is None or self.current[0] != transcript_id:
            return []
        rows = list(self.current[1])
        self.current = next(self.groups, None)
        return rows


def reference_structure(
    transcript: Mapping[str, Any],
    exons: Sequence[DatabaseRow],
    translation: Mapping[str, Any] | None,
) -> TranscriptProteinFeatureResult:
    """Validate a model and derive exon, splice-window and CDS coordinates.

    Genomic positions remain on the reference strand. Pre-mRNA and CDS positions
    follow transcription direction; pre-mRNA includes introns, CDS does not.
    A missing translation produces a noncoding model with no CDS blocks.
    Missing leading codon bases occupy peptide coordinates but have no genome
    coordinates; the first CDS block starts after that offset.
    """
    strand = transcript["strand"]
    start, end = transcript["start"], transcript["end"]
    chromosome, assembly = transcript["chromosome"], transcript["assembly_name"]
    if strand not in {-1, 1} or not 1 <= start <= end or not exons:
        raise ValueError("Invalid transcript strand/span or missing exons")
    # Negative-strand transcripts start at the highest genomic exon coordinate.
    oriented = sorted(exons, key=lambda exon: exon["start"], reverse=strand == -1)
    if [exon["rank"] for exon in oriented] != list(range(1, len(exons) + 1)):
        raise ValueError("Exon ranks disagree with transcript strand/order")
    result: TranscriptProteinFeatureResult = {
        "translation_id": translation["stable_id"] if translation else None,
        "protein_length": None,
        "chromosome": chromosome,
        "strand": strand,
        "transcript_genomic_start": start,
        "transcript_genomic_end": end,
        "premrna_length": end - start + 1,
        "transcript_exons": [],
        "cds_blocks": [],
        "protein_features": [],
        "splice_sites": [],
        "assembly_name": assembly,
    }
    for number, exon in enumerate(oriented, 1):
        lo, hi = exon["start"], exon["end"]
        if not start <= lo <= hi <= end or exon["strand"] != strand:
            raise ValueError("Exon coordinates/strand inconsistent with transcript")
        if exon["seq_region_id"] != transcript["seq_region_id"]:
            raise ValueError("Exon and transcript use different sequence regions")
        if number > 1:
            previous = oriented[number - 2]
            if (strand == 1 and previous["end"] >= lo) or (
                strand == -1 and hi >= previous["start"]
            ):
                raise ValueError("Overlapping exons within transcript")
        premrna_start = lo - start + 1 if strand == 1 else end - hi + 1
        premrna_end = hi - start + 1 if strand == 1 else end - lo + 1
        result["transcript_exons"].append(
            {
                "exon_number": number,
                "chromosome": chromosome,
                "genomic_start": lo,
                "genomic_end": hi,
                "premrna_start": premrna_start,
                "premrna_end": premrna_end,
            }
        )
        for site_type, position, genomic in (
            ("acceptor", premrna_start, lo if strand == 1 else hi),
            ("donor", premrna_end, hi if strand == 1 else lo),
        ):
            # Outer transcript ends are not internal splice junctions. Each retained
            # window spans two exon bases and the two adjacent intron bases.
            if (site_type == "acceptor" and number == 1) or (
                site_type == "donor" and number == len(oriented)
            ):
                continue
            if (strand, site_type) in {(-1, "donor"), (1, "acceptor")}:
                disruption_start, disruption_end = genomic - 2, genomic + 1
            else:
                disruption_start, disruption_end = genomic - 1, genomic + 2
            result["splice_sites"].append(
                {
                    "type": cast(Literal["acceptor", "donor"], site_type),
                    "exon_number": number,
                    "genomic_position": genomic,
                    "premrna_position": position,
                    "disruption_start": disruption_start,
                    "disruption_end": disruption_end,
                }
            )
    if translation is None:
        return result
    exon_by_id = {exon["exon_id"]: exon for exon in oriented}
    first = exon_by_id.get(translation["start_exon_id"])
    last = exon_by_id.get(translation["end_exon_id"])
    if first is None or last is None or first["rank"] > last["rank"]:
        raise ValueError("Translation endpoints do not belong to ordered transcript exons")
    for exon, offset in ((first, translation["seq_start"]), (last, translation["seq_end"])):
        if not 1 <= offset <= exon["end"] - exon["start"] + 1:
            raise ValueError("Translation endpoint offset lies outside exon")
    # Translation offsets are 1-based within the endpoint exons, counted in
    # transcription direction (backwards from exon end on the negative strand).
    first_pos = (
        first["start"] + translation["seq_start"] - 1
        if strand == 1
        else first["end"] - translation["seq_start"] + 1
    )
    last_pos = (
        last["start"] + translation["seq_end"] - 1
        if strand == 1
        else last["end"] - translation["seq_end"] + 1
    )
    if (strand == 1 and first_pos > last_pos) or (strand == -1 and first_pos < last_pos):
        raise ValueError("Translation start follows translation end")
    coding_lo, coding_hi = sorted((first_pos, last_pos))
    # A 5'-incomplete CDS starts partway through the first peptide codon.
    # Ensembl pads that codon with `phase` unknown bases when translating.
    # Keep those peptide coordinates, but map only bases present in the genome:
    # starting at phase + 1 also prevents a missing native start being retained.
    phase = first["phase"]
    if phase not in {-1, 0, 1, 2}:
        raise ValueError("Invalid translation start exon phase")
    start_phase = max(0, phase)
    if start_phase:
        result["cds_start_phase"] = start_phase
    cds_pos = start_phase + 1
    # Clip each exon to the coding span, then concatenate its CDS coordinates.
    for prepared_exon in result["transcript_exons"]:
        lo, hi = (
            max(prepared_exon["genomic_start"], coding_lo),
            min(prepared_exon["genomic_end"], coding_hi),
        )
        if lo > hi:
            continue
        length = hi - lo + 1
        result["cds_blocks"].append(
            {
                "chromosome": chromosome,
                "genomic_start": lo,
                "genomic_end": hi,
                "premrna_start": lo - start + 1 if strand == 1 else end - hi + 1,
                "premrna_end": hi - start + 1 if strand == 1 else end - lo + 1,
                "strand": strand,
                "assembly_name": assembly,
                "cds_start": cds_pos,
                "cds_end": cds_pos + length - 1,
            }
        )
        cds_pos += length
    return result


def protein_segments(
    aa_start: int,
    aa_end: int,
    blocks: Sequence[CDSBlock],
    *,
    block_ends: Sequence[int] | None = None,
) -> list[GenomicSegment]:
    """Map an inclusive amino-acid interval to genomic segments across CDS exons.

    Blocks are ordered by CDS position. Passing their end positions enables a
    binary search; callers mapping many features can reuse that index.
    """
    segments: list[GenomicSegment] = []
    cds_start, cds_end = (aa_start - 1) * 3 + 1, aa_end * 3
    # The preprocessor supplies a once-per-transcript index of ordered CDS blocks.
    if block_ends is not None:
        first = bisect_left(block_ends, cds_start)
        last = bisect_left(block_ends, cds_end) + 1
        selected = blocks[first:last]
    else:
        selected = blocks
    for block in selected:
        lo = max(cds_start, block["cds_start"])
        hi = min(cds_end, block["cds_end"])
        if lo > hi:
            continue
        start_offset = lo - block["cds_start"]
        end_offset = hi - block["cds_start"]
        # Return ascending genomic bounds even when CDS direction is reversed.
        if block["strand"] == 1:
            start = block["genomic_start"] + start_offset
            end = block["genomic_start"] + end_offset
        else:
            start = block["genomic_end"] - end_offset
            end = block["genomic_end"] - start_offset
        segments.append(
            {
                "chromosome": block["chromosome"],
                "start": start,
                "end": end,
                "strand": block["strand"],
                "assembly_name": block["assembly_name"],
            }
        )
    return segments


# Source compatibility checks and bulk protein annotation metadata.


def protein_analyses(db: sqlite3.Connection) -> AnalysisMetadata:
    """Use database names for feature sources, not Ensembl pipeline names.

    For example, logic_name='hmmpanther' describes how the analysis ran;
    db='PANTHER' identifies the annotation source stored with each feature.
    Small offline fixtures may have only logic_name, so retain that fallback.
    Structure pipelines override database labels so their mappings can always
    be excluded from functional annotations.
    """
    analyses = {}
    for row in dict_rows(db, "SELECT * FROM ensembl_analysis"):
        source = row.get("db") or row["logic_name"]
        logic_name = (row["logic_name"] or "").lower()
        if logic_name in STRUCTURE_SOURCES:
            source = logic_name
        if (source or "").lower() in {"panther", "hmmpanther"}:
            source = "PANTHER"
        analyses[row["analysis_id"]] = (source, row.get("db_version"))
    return analyses


def download_panther_classifications(version: str, cache: Path) -> tuple[Path, str, str]:
    """Fetch the exact release, accepting the historical trailing underscore.

    Prefer an already cached filename so archived files also work offline.
    Only a 404 permits the alternate name; other HTTP/network errors propagate.
    """
    directory = cache / version
    root = PANTHER_CLASSIFICATIONS_URL + version + "/PANTHER_Sequence_Classification_files/"
    names = [f"PTHR{version}_human", f"PTHR{version}_human_"]
    names.sort(
        key=lambda name: (
            not ((directory / name).is_file() and (directory / (name + ".sha256")).is_file())
        )
    )
    for index, name in enumerate(names):
        url = root + name
        try:
            path, digest = download(url, directory)
            return path, digest, url
        except urllib.error.HTTPError as exc:
            if exc.code != 404 or index == len(names) - 1:
                raise
            LOG.info(
                "PANTHER %s filename unavailable; trying the alternate archived filename", version
            )
    raise AssertionError("Unreachable")


def panther_release(analyses: AnalysisMetadata) -> str | None:
    """Require one exact PANTHER release when automatic metadata is needed."""
    versions = {version for source, version in analyses.values() if source == "PANTHER"}
    if not versions:
        return None
    if len(versions) != 1 or not re.fullmatch(r"\d+(?:\.\d+)*", next(iter(versions)) or ""):
        raise ValueError(
            "Cannot determine one PANTHER release from Ensembl analysis metadata; "
            "supply --panther-classifications FILE"
        )
    return next(iter(versions))


def panther_entry_rows(path: Path) -> Iterator[tuple[str, str, str]]:
    """Read human UniProt-to-subfamily assignments from PANTHER's bulk TSV.

    Older files place the subfamily ID/name in columns 3/5; newer files use
    columns 4/6. Detect the layout from the classification ID, independently of
    the filename. Protein accessions always come from the identifier in column 1.
    Family/subfamily classifications do not supply new domain coordinates.
    """
    with input_progress(path, "Import PANTHER", compressed=False, text=True) as (
        stream,
        bar,
        update,
    ):
        count = 0
        columns = None
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5 or not fields[0].startswith("HUMAN|"):
                raise ValueError(f"Invalid human PANTHER classification at {path}:{number}")
            accession = re.search(r"(?:^|\|)UniProtKB[:=]([^|]+)", fields[0])
            if accession is None:
                raise ValueError(f"Missing UniProt accession at {path}:{number}")
            if columns is None:
                layouts = [
                    (id_col, name_col)
                    for id_col, name_col in ((3, 5), (2, 4))
                    if len(fields) > name_col
                ]
                columns = next(
                    (
                        layout
                        for layout in layouts
                        if re.fullmatch(r"PTHR\d+:SF\d+", fields[layout[0]])
                    ),
                    None,
                )
                # Unassigned leading records can still identify an empty ID
                # column. Once detected, keep the layout fixed to catch bad rows.
                if columns is None:
                    columns = next((layout for layout in layouts if not fields[layout[0]]), None)
                if columns is None:
                    raise ValueError(
                        f"Invalid PANTHER subfamily at {path}:{number}; "
                        "expected an ID in column 3 or 4"
                    )
            id_col, name_col = columns
            if len(fields) <= name_col:
                raise ValueError(f"Invalid human PANTHER classification at {path}:{number}")
            subfamily, name = fields[id_col], fields[name_col]
            # Records without a named subfamily cannot support the family fallback.
            if not subfamily:
                continue
            if not re.fullmatch(r"PTHR\d+:SF\d+", subfamily):
                raise ValueError(f"Invalid PANTHER subfamily at {path}:{number}: {subfamily}")
            if not name:
                continue
            count += 1
            yield accession[1], subfamily, name
            if number % 1000 == 0:
                bar.set_postfix(records=count, refresh=False)
                update()
        if not count:
            raise ValueError(f"No named human PANTHER subfamilies in {path}")
        bar.set_postfix(records=count, refresh=False)


def panther_translation_annotations(db: sqlite3.Connection) -> dict[int, tuple[str, str]]:
    """Match exact protein cross-references; never propagate by gene name.

    An ambiguous protein assignment is omitted. The caller also requires its
    family to agree with the coordinate-bearing Ensembl PANTHER hit.
    """
    tables = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "ensembl_object_xref" not in tables:
        return {}
    matches: dict[int, set[tuple[str, str]]] = {}
    # Start from the small human lookup and use the imported accession/xref
    # indexes, rather than scanning every row in the much larger core tables.
    for translation_id, subfamily, name in db.execute("""
            SELECT o.ensembl_id, p.subfamily_id, p.name
            FROM ff_panther p CROSS JOIN ensembl_xref x ON x.dbprimary_acc=p.accession
            CROSS JOIN ensembl_object_xref o ON o.xref_id=x.xref_id
            WHERE o.ensembl_object_type='Translation'
        """):
        matches.setdefault(translation_id, set()).add((subfamily, name))
    return {key: next(iter(values)) for key, values in matches.items() if len(values) == 1}


def validate_reference_source(db: sqlite3.Connection) -> dict[str, str]:
    """Reject incompatible or incomplete source databases before changing them."""
    available = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    required = {"ensembl_" + name for name in PREPROCESS_TABLES}
    required.update({"build_metadata", "sequences", "sequence_chunks"})
    missing = required - available
    if missing:
        if (
            "build_metadata" in available
            and dict(db.execute("SELECT key, value FROM build_metadata")).get("reference_kind")
            == "runtime"
        ):
            raise ValueError(
                "This compact prebuilt reference has no Ensembl source tables. "
                "Install a newer prebuilt with prepare-data --release RELEASE --force, "
                "or build a full reference with --from-source before reprocessing."
            )
        raise ValueError(
            f"Preprocessing needs {sorted(missing)}; rebuild with 'fusion-function prepare-data --from-source'"
        )
    metadata = dict(db.execute("SELECT key, value FROM build_metadata"))
    if metadata.get("format_version") != "1":
        raise ValueError(
            "Unsupported reference format; rebuild with 'fusion-function prepare-data'"
        )
    if metadata.get("species") != SPECIES:
        raise ValueError("Only human reference databases are supported")
    if (
        metadata.get("assembly") != ASSEMBLY
        or not db.execute(
            "SELECT 1 FROM ensembl_coord_system WHERE version=? LIMIT 1", (ASSEMBLY,)
        ).fetchone()
    ):
        raise ValueError("Only GRCh38 reference databases are supported")
    if not re.fullmatch(r"[1-9]\d*", metadata.get("release", "")):
        raise ValueError("Reference metadata must contain a positive Ensembl release number")
    if (
        metadata.get("sequence_chunk_size") != str(CHUNK_SIZE)
        or metadata.get("sequence_codec") != "zlib"
    ):
        raise ValueError("Unsupported sequence chunk format")
    if metadata.get("preprocessing_version") not in {None, "1"}:
        raise ValueError("Unsupported preprocessing version")
    kinds = {row[0] for row in db.execute("SELECT DISTINCT kind FROM sequences")}
    if not set(SEQUENCE_TYPES) <= kinds:
        raise ValueError(
            "Preprocessing requires DNA and peptide sequences; rebuild with 'fusion-function prepare-data'"
        )
    return metadata


def interpro_entry_rows(path: Path) -> Iterator[tuple[str, str | None, str | None]]:
    """Read current and archived InterPro entry lists with the same validation."""
    with input_progress(path, "Import InterPro", compressed=False, text=True) as (
        stream,
        bar,
        update,
    ):
        header = stream.readline().rstrip("\r\n").split("\t")
        expected = ["ENTRY_AC", "ENTRY_TYPE", "ENTRY_NAME"]
        if not set(expected) <= set(header):
            raise ValueError("InterPro TSV must have ENTRY_AC, ENTRY_TYPE and ENTRY_NAME columns")
        indices = [header.index(name) for name in expected]
        for number, line in enumerate(stream, 2):
            if not line.strip():
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                raise ValueError(f"InterPro TSV line {number}: wrong field count")
            accession, entry_type, name = [values[index] for index in indices]
            if not re.fullmatch(r"IPR\d+", accession):
                raise ValueError(f"Invalid InterPro accession on line {number}")
            entry_type = entry_type.strip().lower().replace(" ", "_").replace("-", "_")
            yield accession, name or None, entry_type or None
            if number % 1000 == 0:
                update()


def restore_interpro_history(
    db: sqlite3.Connection, metadata: InterProMetadata, cache: Path
) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Supplement missing types, newest archive first; never overwrite current types.

    Explicit local entry lists do not trigger archive downloads.
    """
    # This small mapping table includes all possible core accessions. Only if
    # archives cannot resolve some do we scan the much larger feature table.
    required = {row[0] for row in db.execute("SELECT DISTINCT interpro_ac FROM ensembl_interpro")}
    missing = {accession for accession in required if not metadata.get(accession, (None, None))[1]}
    if not missing:
        return {}, []
    LOG.warning(
        "%s Ensembl InterPro accessions lack current metadata; searching archived entry lists before transcript preprocessing",
        len(missing),
    )
    versions = sorted(
        (name for name in listing(INTERPRO_RELEASES_URL) if re.fullmatch(r"\d+(?:\.\d+)?", name)),
        key=lambda name: tuple(int(part) for part in name.split(".")),
        reverse=True,
    )
    restored, sources = {}, []
    for version in versions:
        url = INTERPRO_RELEASES_URL + version + "/entry.list"
        try:
            path, digest = download(url, cache / "releases" / version)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                continue
            raise
        recovered = []
        for accession, name, entry_type in interpro_entry_rows(path):
            if accession in missing and entry_type:
                recovered.append((accession, name, entry_type))
                metadata[accession] = name or metadata.get(accession, (None, None))[0], entry_type
                restored[accession] = version
                missing.remove(accession)
        if recovered:
            db.executemany(
                "INSERT INTO ff_interpro VALUES (?, ?, ?) ON CONFLICT(interpro_id) "
                "DO UPDATE SET name=COALESCE(excluded.name, ff_interpro.name), "
                "entry_type=excluded.entry_type",
                recovered,
            )
            sources.append({"release": version, "url": url, "sha256": digest})
            LOG.info(
                "Recovered %s historical InterPro entries from release %s; %s remaining",
                len(recovered),
                version,
                len(missing),
            )
        if not missing:
            break
    if missing:
        analyses = protein_analyses(db)
        used = {
            row[0]
            for row in sqlite_phase(
                db,
                """
            SELECT DISTINCT i.interpro_ac, pf.analysis_id FROM ensembl_protein_feature pf
            JOIN ensembl_translation tr USING (translation_id)
            JOIN ensembl_transcript t ON t.transcript_id=tr.transcript_id
                                     AND t.canonical_translation_id=tr.translation_id
            JOIN ensembl_interpro i ON i.id=pf.hit_name
            WHERE t.is_current=1
        """,
                "Check unresolved InterPro metadata",
            )
            if (analyses.get(row[1], (None, None))[0] or "").lower() not in STRUCTURE_SOURCES
        }
        missing.intersection_update(used)
    if missing:
        examples = ", ".join(sorted(missing)[:10])
        raise ValueError(
            f"Missing InterPro entry types for {len(missing)} accession(s) after checking archived metadata: "
            f"{examples}. Supply a complete --interpro-entries FILE; transcript preprocessing has not started."
        )
    return restored, sources


# Derived-table preparation: one transaction, one transcript at a time.


def _record_transcript_error(
    summary: dict[str, TranscriptErrorSummary], transcript_id: str, message: str
) -> None:
    """Group errors by reason while keeping variable details in the stored payload."""
    # CDS/peptide mismatch lengths and unsupported assembly names vary, so group
    # by the fixed message prefix. Each original full error remains retrievable.
    reason = message.partition(":")[0]
    entry = summary.setdefault(reason, {"count": 0, "example_transcripts": []})
    entry["count"] += 1
    if len(entry["example_transcripts"]) < 3:
        entry["example_transcripts"].append(transcript_id)


def uniprot_lookup_fingerprint(
    db: sqlite3.Connection, metadata: Mapping[str, str], input_sha256: str
) -> str | None:
    """Key verified lookups by XML, mapping code and immutable imported sources.

    Source checksums describe the raw tables and sequences imported by the
    builder. References without this provenance are reverified, never reused.
    """
    if not db.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='source_files'"
    ).fetchone():
        return None
    sources = list(db.execute("SELECT url, sha256, bytes FROM source_files ORDER BY url"))
    if not sources or not all(row[1] for row in sources):
        return None
    fingerprint = {
        "uniprot_sha256": input_sha256,
        "mapping_sha256": sha256_file(Path(__file__).with_name("uniprot.py")),
        "ensembl_sources": sources,
        "reference": {
            key: metadata.get(key)
            for key in (
                "format_version",
                "release",
                "species",
                "assembly",
                "sequence_codec",
                "sequence_chunk_size",
            )
        },
    }
    return hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()


def encode_transcript_payload(payload: ProteinFeatureResponse) -> str | bytes:
    """Compress large models while keeping small error records SQL-readable.

    Level 1 limits preprocessing CPU cost. This is lossless storage compression;
    it does not remove features or change the annotation API's returned model.
    """
    encoded = json.dumps(payload, separators=(",", ":"))
    return encoded if "error" in payload else zlib.compress(encoded.encode("utf-8"), level=1)


def decode_transcript_payload(encoded: str | bytes) -> ProteinFeatureResponse:
    """Decode compressed models or plain JSON errors and existing references."""
    if isinstance(encoded, bytes):
        encoded = zlib.decompress(encoded).decode("utf-8")
    return cast("ProteinFeatureResponse", json.loads(encoded))


def preprocess_reference(
    db: sqlite3.Connection,
    interpro_entries: Path | None = None,
    *,
    interpro_archive_cache: Path | None = None,
    panther_classifications: Path | None = None,
    panther_cache: Path | None = None,
    uniprot_features: Path | None = None,
    uniprot_source: dict[str, str] | None = None,
) -> dict[str, int]:
    """Build all derived tables in one transaction; leave prior tables on failure.

    Requires raw annotation tables and DNA/peptide FASTA. Invalid individual
    transcript models are recorded as errors, rather than silently corrected.
    No nucleotide sequences are duplicated in the derived transcript records.
    """
    from .uniprot import CuratedFeature, import_features

    metadata = validate_reference_source(db)
    analyses = protein_analyses(db)
    panther_source = None
    if panther_classifications is None and panther_cache is not None:
        version = panther_release(analyses)
        if version:
            panther_classifications, digest, url = download_panther_classifications(
                version, panther_cache
            )
            panther_source = {"version": version, "url": url, "sha256": digest}
    if uniprot_features is not None and uniprot_source is None:
        uniprot_source = {
            "path": str(uniprot_features.resolve()),
            "sha256": sha256_file(uniprot_features),
        }
    uniprot_fingerprint = (
        uniprot_lookup_fingerprint(db, metadata, uniprot_source["sha256"])
        if uniprot_features is not None and uniprot_source is not None
        else None
    )
    if panther_classifications is not None and panther_source is None:
        panther_source = {
            "path": str(panther_classifications.resolve()),
            "sha256": sha256_file(panther_classifications),
        }
    LOG.info("Preprocessing reference transcript models, domains and splice sites")
    db.commit()
    # These are reproducible public reference annotations. Some SQLite builds
    # default to zeroing deleted content, which rewrites entire large tables and
    # their rollback journals. Reclaim pages without scrubbing; journaling and
    # atomic rollback remain enabled. Restore the caller's setting on every exit.
    secure_delete = int(db.execute("PRAGMA main.secure_delete").fetchone()[0])
    # FAST reads back as 2, but setting the numeric value 2 means ON. Restore
    # the keyword so a caller using FAST keeps that exact policy.
    previous_secure_delete = ("OFF", "ON", "FAST")[secure_delete]
    db.execute("PRAGMA main.secure_delete=OFF")
    # Derived-table replacement, including DDL, is atomic. A failed metadata
    # lookup or model build rolls back to the prior usable reference tables.
    transcript_progress = None
    try:
        LOG.info(
            "SQLite reference preprocessing: secure_delete=%s (previously %s), "
            "journal_mode=%s, synchronous=%s, auto_vacuum=%s, page_size=%s",
            db.execute("PRAGMA main.secure_delete").fetchone()[0],
            previous_secure_delete,
            db.execute("PRAGMA main.journal_mode").fetchone()[0],
            db.execute("PRAGMA main.synchronous").fetchone()[0],
            db.execute("PRAGMA main.auto_vacuum").fetchone()[0],
            db.execute("PRAGMA main.page_size").fetchone()[0],
        )
        db.execute("BEGIN IMMEDIATE")
        if uniprot_features is not None:
            if (
                uniprot_fingerprint is not None
                and uniprot_fingerprint == metadata.get("uniprot_lookup_fingerprint")
                and metadata.get("uniprot_counts") is not None
                and db.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='ff_uniprot_features'"
                ).fetchone()
            ):
                LOG.info("Reusing sequence-verified UniProt lookups: inputs and mapping unchanged")
                uniprot_counts = json.loads(metadata["uniprot_counts"])
            else:
                uniprot_counts = import_features(db, uniprot_features)
        else:
            uniprot_counts = None
        # An offline low-level call can reuse already verified UniProt lookups.
        db.execute(
            "CREATE TABLE IF NOT EXISTS ff_uniprot_features (sequence_id TEXT NOT NULL, "
            "feature_id TEXT NOT NULL, payload TEXT NOT NULL, "
            "PRIMARY KEY(sequence_id, feature_id)) WITHOUT ROWID"
        )
        # Persist the bulk lookup so subsequent offline preprocessing can reuse it.
        # Import and all derived changes participate in the same rollback boundary.
        db.execute(
            "CREATE TABLE IF NOT EXISTS ff_panther (accession TEXT PRIMARY KEY, "
            "subfamily_id TEXT NOT NULL, name TEXT NOT NULL) WITHOUT ROWID"
        )
        if panther_classifications is not None:
            db.execute("DELETE FROM ff_panther")
            db.executemany(
                "INSERT INTO ff_panther VALUES (?, ?, ?)",
                panther_entry_rows(panther_classifications),
            )
        with sqlite_activity("Match PANTHER subfamilies to Ensembl proteins"):
            panther_annotations = (
                panther_translation_annotations(db)
                if db.execute("SELECT 1 FROM ff_panther LIMIT 1").fetchone()
                else {}
            )
        LOG.info(
            "PANTHER subfamily assignments matched %s Ensembl translations by protein cross-reference",
            f"{len(panther_annotations):,}",
        )
        with sqlite_activity("Preserve existing InterPro metadata"):
            has_interpro = db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='ff_interpro'"
            ).fetchone()
            previous_interpro = (
                list(db.execute("SELECT interpro_id, name, entry_type FROM ff_interpro"))
                if has_interpro
                else []
            )
        for table in ("ff_transcripts", "ff_interpro"):
            with sqlite_activity("Drop previous " + table):
                db.execute("DROP TABLE IF EXISTS " + sql_name(table))
        with sqlite_activity("Create derived transcript and InterPro tables"):
            # Large model payloads belong in ordinary rowid-table leaves, not
            # the intermediate nodes of a WITHOUT ROWID primary-key tree.
            # Internal transcript-ID order then appends model rows sequentially;
            # only the much smaller stable-ID index needs random inserts.
            db.execute(
                "CREATE TABLE ff_transcripts (transcript_id TEXT NOT NULL PRIMARY KEY, version INTEGER, "
                "translation_id TEXT, translation_version INTEGER, status TEXT NOT NULL, "
                "payload BLOB NOT NULL)"
            )
            db.execute(
                "CREATE TABLE ff_interpro (interpro_id TEXT PRIMARY KEY, name TEXT, entry_type TEXT) WITHOUT ROWID"
            )
        # Priority: core names < previously stored types < supplied/current list.
        # History fills only remaining gaps; it never overwrites a known type.
        with sqlite_activity("Load InterPro names from Ensembl cross-references"):
            db.execute(
                "INSERT INTO ff_interpro SELECT dbprimary_acc, "
                "COALESCE(MAX(NULLIF(description,'')), MAX(NULLIF(display_label,''))), NULL "
                "FROM ensembl_xref WHERE dbprimary_acc LIKE 'IPR%' GROUP BY dbprimary_acc"
            )
            db.executemany(
                "INSERT INTO ff_interpro VALUES (?, ?, ?) ON CONFLICT(interpro_id) "
                "DO UPDATE SET name=COALESCE(excluded.name, ff_interpro.name), "
                "entry_type=excluded.entry_type",
                previous_interpro,
            )
        provided_accessions = set()
        if interpro_entries:
            LOG.info("Importing InterPro entry names/types from %s", interpro_entries)
            for row in interpro_entry_rows(interpro_entries):
                provided_accessions.add(row[0])
                db.execute(
                    "INSERT INTO ff_interpro VALUES (?, ?, ?) ON CONFLICT(interpro_id) "
                    "DO UPDATE SET name=excluded.name, entry_type=excluded.entry_type",
                    row,
                )
        ipr_metadata = {row[0]: (row[1], row[2]) for row in db.execute("SELECT * FROM ff_interpro")}
        historical_entries = json.loads(metadata.get("interpro_historical_entries", "{}"))
        historical_sources = json.loads(metadata.get("interpro_historical_sources", "[]"))
        # Entries explicitly present in the new list take precedence over history.
        for accession in provided_accessions:
            historical_entries.pop(accession, None)
        if interpro_archive_cache is not None:
            restored, sources = restore_interpro_history(db, ipr_metadata, interpro_archive_cache)
            historical_entries.update(restored)
            historical_sources.extend(
                source for source in sources if source not in historical_sources
            )
        with sqlite_activity("Read genome and protein sequence lengths"):
            genome_lengths = dict(
                db.execute("SELECT sequence_id, length FROM sequences WHERE kind='dna'")
            )
            protein_lengths: dict[str, list[tuple[str, int]]] = {}
            for stable_id, sequence_id, length in db.execute(
                "SELECT stable_id, sequence_id, length FROM sequences WHERE kind='pep'"
            ):
                protein_lengths.setdefault(stable_id, []).append((sequence_id, length))
        curated_features = GroupStream(
            dict_rows(
                db,
                """
            SELECT t.transcript_id, u.payload
            FROM ensembl_transcript t JOIN ensembl_translation tr
              ON tr.translation_id=t.canonical_translation_id
            JOIN ff_uniprot_features u ON u.sequence_id IN
              (tr.stable_id, tr.stable_id || '.' || tr.version)
            WHERE t.is_current=1 ORDER BY t.transcript_id, u.feature_id
        """,
                description="Prepare ordered UniProt feature stream",
            )
        )
        # All three queries use the same ascending internal ID order. GroupStream
        # avoids per-transcript SQL queries and loading all exon/features into RAM.
        exons = GroupStream(
            dict_rows(
                db,
                """
            SELECT et.transcript_id, et.rank, e.exon_id, e.seq_region_id,
                   e.seq_region_start AS start, e.seq_region_end AS end,
                   e.seq_region_strand AS strand, e.phase
            FROM ensembl_exon_transcript et JOIN ensembl_exon e USING (exon_id)
            ORDER BY et.transcript_id, et.rank
        """,
                description="Prepare ordered exon stream",
            )
        )
        # Restrict hits to each current transcript's canonical translation; using
        # alternate translations here would mix incompatible peptide coordinates.
        features = GroupStream(
            dict_rows(
                db,
                """
            SELECT tr.transcript_id, pf.protein_feature_id, pf.seq_start AS start,
                   pf.seq_end AS end, pf.hit_name AS feature_id,
                   pf.hit_description AS description, a.logic_name AS source, pf.analysis_id,
                   i.interpro_ac AS interpro_id
            FROM ensembl_protein_feature pf
            JOIN ensembl_translation tr USING (translation_id)
            JOIN ensembl_transcript t ON t.transcript_id=tr.transcript_id
                                     AND t.canonical_translation_id=tr.translation_id
            LEFT JOIN ensembl_analysis a USING (analysis_id)
            LEFT JOIN ensembl_interpro i ON i.id=pf.hit_name
            WHERE t.is_current=1
            ORDER BY tr.transcript_id, pf.protein_feature_id, i.interpro_ac
        """,
                description="Prepare ordered protein feature stream",
            )
        )
        transcripts = dict_rows(
            db,
            """
            SELECT t.transcript_id AS internal_id, t.stable_id AS transcript_id,
                   t.version, t.seq_region_id, t.seq_region_start AS start,
                   t.seq_region_end AS end, t.seq_region_strand AS strand,
                   t.canonical_translation_id, r.name AS chromosome,
                   cs.version AS assembly_name, tr.stable_id AS translation_id,
                   tr.version AS translation_version, tr.seq_start, tr.start_exon_id,
                   tr.seq_end, tr.end_exon_id,
                   EXISTS(SELECT 1 FROM ensembl_translation alt
                          WHERE alt.transcript_id=t.transcript_id) AS has_translation
            FROM ensembl_transcript t JOIN ensembl_seq_region r USING (seq_region_id)
            JOIN ensembl_coord_system cs USING (coord_system_id)
            LEFT JOIN ensembl_translation tr ON tr.translation_id=t.canonical_translation_id
            WHERE t.is_current=1 ORDER BY t.transcript_id
        """,
            description="Prepare transcript stream",
        )
        counts = {"ready": 0, "noncoding": 0, "error": 0, "protein_features": 0}
        errors: dict[str, TranscriptErrorSummary] = {}
        missing_entry_types: set[str] = set()
        total = sqlite_phase(
            db,
            "SELECT COUNT(*) FROM ensembl_transcript t "
            "JOIN ensembl_seq_region r USING (seq_region_id) "
            "JOIN ensembl_coord_system cs USING (coord_system_id) "
            "WHERE t.is_current=1",
            "Count current transcripts",
        )[0][0]
        transcript_progress = progress(
            transcripts, total=total, desc="Preprocess transcripts", unit="transcript"
        )
        payload: ProteinFeatureResponse
        for number, transcript in enumerate(transcript_progress, 1):
            exon_rows = exons.take(transcript["internal_id"])
            feature_rows = features.take(transcript["internal_id"])
            curated_rows = curated_features.take(transcript["internal_id"])
            translation = None
            if transcript["canonical_translation_id"] is not None:
                translation = {**transcript, "stable_id": transcript["translation_id"]}
            try:
                if transcript["assembly_name"] != ASSEMBLY:
                    raise ValueError(
                        f"Unsupported transcript coordinate system: {transcript['assembly_name']}"
                    )
                if translation is not None and not translation["stable_id"]:
                    raise ValueError("Canonical translation is missing from source table")
                if translation is None and transcript["has_translation"]:
                    raise ValueError("Transcript has translations but no canonical translation")
                payload = reference_structure(transcript, exon_rows, translation)
                if transcript["chromosome"] not in genome_lengths:
                    raise ValueError("Transcript sequence region is absent from genomic FASTA")
                if transcript["end"] > genome_lengths[transcript["chromosome"]]:
                    raise ValueError("Transcript span exceeds genomic FASTA sequence")
                status = "noncoding" if translation is None else "ready"
                if translation:
                    candidates = protein_lengths.get(translation["stable_id"], [])
                    expected_id = f"{translation['stable_id']}.{transcript['translation_version']}"
                    candidates = [
                        record
                        for record in candidates
                        if record[0] in {expected_id, translation["stable_id"]}
                    ]
                    if len(candidates) != 1:
                        raise ValueError(
                            "Canonical translation peptide is missing/ambiguous or has a different version"
                        )
                    protein_length = candidates[0][1]
                    payload["protein_length"] = protein_length
                    cds_length = sum(
                        block["cds_end"] - block["cds_start"] + 1 for block in payload["cds_blocks"]
                    )
                    # Include Ensembl's virtual leading codon bases in the length
                    # comparison, without adding them to genomic CDS blocks.
                    # Peptide FASTA omits the optional terminal stop codon.
                    start_phase = payload.get("cds_start_phase", 0)
                    if cds_length + start_phase not in {protein_length * 3, protein_length * 3 + 3}:
                        raise ValueError(
                            f"CDS/peptide length mismatch: {cds_length} bp, {protein_length} aa, "
                            f"start phase {start_phase}"
                        )
                    block_ends = [block["cds_end"] for block in payload["cds_blocks"]]
                    # Multiple signatures/InterPro mappings can share an interval.
                    # Reuse coordinates within this transcript, then discard them.
                    segment_cache: dict[tuple[int, int], tuple[list[GenomicSegment], int, int]] = {}
                    for feature in feature_rows:
                        source = analyses.get(feature["analysis_id"], (feature["source"], None))[0]
                        if (source or "").lower() in STRUCTURE_SOURCES:
                            # Whole-protein structure mappings are not domains/sites.
                            # Exclude them before validating functional coordinates:
                            # their intervals may exceed this transcript's peptide.
                            continue
                        if not 1 <= feature["start"] <= feature["end"] <= protein_length:
                            raise ValueError(
                                f"Protein feature falls outside peptide: {source} "
                                f"{feature['feature_id']} {feature['start']}-{feature['end']}, "
                                f"peptide length {protein_length}"
                            )
                        interval = feature["start"], feature["end"]
                        mapped = segment_cache.get(interval)
                        if mapped is None:
                            segments = protein_segments(
                                *interval, payload["cds_blocks"], block_ends=block_ends
                            )
                            mapped = (
                                segments,
                                min(segment["start"] for segment in segments),
                                max(segment["end"] for segment in segments),
                            )
                            segment_cache[interval] = mapped
                        segments, genomic_start, genomic_end = mapped
                        name, entry_type = ipr_metadata.get(feature["interpro_id"], (None, None))
                        subfamily = (
                            feature["feature_id"]
                            if (source or "").upper() == "PANTHER"
                            and ":SF" in (feature["feature_id"] or "")
                            else None
                        )
                        subfamily_name = feature["description"] if subfamily else None
                        if source == "PANTHER":
                            assignment = panther_annotations.get(
                                transcript["canonical_translation_id"]
                            )
                            hit = feature["feature_id"] or ""
                            if assignment and assignment[0].split(":")[0] == hit.split(":")[0]:
                                if subfamily is None or subfamily == assignment[0]:
                                    subfamily, subfamily_name = assignment
                        payload["protein_features"].append(
                            {
                                "feature_id": feature["feature_id"] or None,
                                "source": source or None,
                                "interpro_id": feature["interpro_id"],
                                "start": feature["start"],
                                "end": feature["end"],
                                "cds_start": (feature["start"] - 1) * 3 + 1,
                                "cds_end": feature["end"] * 3,
                                "chromosome": transcript["chromosome"],
                                "strand": transcript["strand"],
                                "assembly_name": transcript["assembly_name"],
                                "genomic_start": genomic_start,
                                "genomic_end": genomic_end,
                                "genomic_segments": segments,
                                "description": feature["description"] or name,
                                "interpro_name": name,
                                "interpro_entry_type": entry_type,
                                "panther_subfamily_id": subfamily,
                                "panther_subfamily_description": subfamily_name,
                            }
                        )
                    for curated_row in curated_rows:
                        curated = cast(CuratedFeature, json.loads(curated_row["payload"]))
                        start, end = curated["start"], curated["end"]
                        if not 1 <= start <= end <= protein_length:
                            raise ValueError("Verified UniProt feature falls outside peptide")
                        segments = protein_segments(
                            start, end, payload["cds_blocks"], block_ends=block_ends
                        )
                        payload["protein_features"].append(
                            {
                                "feature_id": curated["feature_id"],
                                "feature_type": curated["feature_type"],
                                "description": curated["description"],
                                "start": start,
                                "end": end,
                                "uniprot_accession": curated["uniprot_accession"],
                                "uniprot_isoform": curated["uniprot_isoform"],
                                "evidence": curated["evidence"],
                                "source": "UniProtKB",
                                "interpro_id": None,
                                "cds_start": (start - 1) * 3 + 1,
                                "cds_end": end * 3,
                                "chromosome": transcript["chromosome"],
                                "strand": transcript["strand"],
                                "assembly_name": transcript["assembly_name"],
                                "genomic_start": min(segment["start"] for segment in segments),
                                "genomic_end": max(segment["end"] for segment in segments),
                                "genomic_segments": segments,
                                "interpro_name": None,
                                "interpro_entry_type": None,
                                "panther_subfamily_id": None,
                                "panther_subfamily_description": None,
                            }
                        )
                    counts["protein_features"] += len(payload["protein_features"])
            except ValueError as exc:
                # Preserve the reason for unsupported models so readers can explain
                # the failure, rather than reporting that the transcript is absent.
                status, payload = "error", {"error": f"{transcript['transcript_id']}: {exc}"}
                _record_transcript_error(errors, transcript["transcript_id"], str(exc))
            if status == "ready":
                payload = cast("TranscriptProteinFeatureResult", payload)
                missing_entry_types.update(
                    feature["interpro_id"]
                    for feature in payload["protein_features"]
                    if feature["interpro_id"] and not feature["interpro_entry_type"]
                )
            counts[status] += 1
            db.execute(
                "INSERT INTO ff_transcripts VALUES (?, ?, ?, ?, ?, ?)",
                (
                    transcript["transcript_id"],
                    transcript["version"],
                    transcript["translation_id"],
                    transcript["translation_version"],
                    status,
                    encode_transcript_payload(payload),
                ),
            )
            if number % 1000 == 0:
                transcript_progress.set_postfix(
                    ready=counts["ready"], errors=counts["error"], refresh=False
                )
        if missing_entry_types:
            examples = ", ".join(sorted(missing_entry_types)[:10])
            raise ValueError(
                f"Missing InterPro entry types for {len(missing_entry_types)} accession(s): {examples}. "
                "Supply a complete --interpro-entries FILE or fetch compatible metadata; "
                "the previous reference has been preserved."
            )
        with sqlite_activity("Index prepared transcripts"):
            db.execute("CREATE INDEX ff_transcript_translation ON ff_transcripts (translation_id)")
            db.execute("CREATE INDEX ff_transcript_status ON ff_transcripts (status)")
        # Keep provenance in the same transaction as the tables it describes.
        errors = dict(sorted(errors.items(), key=lambda item: (-item[1]["count"], item[0])))
        preprocessing_metadata = {
            "preprocessing_version": "1",
            "feature_annotation_version": "3",
            "cds_mapping_version": "2",
            "transcript_payload_codec": TRANSCRIPT_PAYLOAD_CODEC,
            "preprocessing_counts": json.dumps(counts),
            "preprocessing_errors": json.dumps(errors),
            "preprocessed_utc": datetime.now(timezone.utc).isoformat(),
            "interpro_entries_sha256": (
                sha256_file(interpro_entries)
                if interpro_entries
                else metadata.get("interpro_entries_sha256", "")
            ),
            "interpro_entry_types_loaded": str(
                any(value[1] for value in ipr_metadata.values())
            ).lower(),
            "interpro_historical_entries": json.dumps(historical_entries, sort_keys=True),
            "interpro_historical_sources": json.dumps(historical_sources, sort_keys=True),
            "interpro_preserved_metadata_sha256": (
                metadata.get("interpro_entries_sha256", "")
                if interpro_entries and previous_interpro
                else metadata.get("interpro_preserved_metadata_sha256", "")
            ),
        }
        if uniprot_source is not None:
            preprocessing_metadata["uniprot_features_source"] = json.dumps(
                uniprot_source, sort_keys=True
            )
            preprocessing_metadata["uniprot_counts"] = json.dumps(uniprot_counts, sort_keys=True)
            # An import without recorded source provenance invalidates any
            # prior reuse key. Empty values never permit reuse.
            preprocessing_metadata["uniprot_lookup_fingerprint"] = uniprot_fingerprint or ""
        if panther_source is not None:
            preprocessing_metadata["panther_classifications_source"] = json.dumps(
                panther_source, sort_keys=True
            )
        for key, value in preprocessing_metadata.items():
            db.execute(
                "INSERT INTO build_metadata VALUES (?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (key, value),
            )
        if counts["ready"] == 0:
            raise ValueError(
                "Preprocessing produced no usable coding transcripts; check source tables and FASTA"
            )
        with sqlite_activity("Commit preprocessed reference"):
            db.commit()
    except BaseException:
        with sqlite_activity("Roll back preprocessing transaction"):
            db.rollback()
        raise
    finally:
        try:
            if transcript_progress is not None:
                transcript_progress.close()
        finally:
            db.execute(f"PRAGMA main.secure_delete={previous_secure_delete}")
    LOG.info("Preprocessing complete: %s", counts)
    if counts["error"]:
        LOG.warning(
            "%s transcripts have stored errors; inspect ff_transcripts WHERE status='error'",
            counts["error"],
        )
        for reason, summary in errors.items():
            LOG.warning(
                "  %s: %s transcript(s); examples: %s",
                reason,
                f"{summary['count']:,}",
                ", ".join(summary["example_transcripts"]),
            )
    return counts


# Read-only runtime access to a completed reference.


class ReferenceReader:
    """Read-only local reference access; reuse one reader across fusion calls.

    get_transcript() returns prepared coordinates, splice sites and features,
    reconstructing pre-mRNA from the genome when requested.
    Genome chunk caching is LRU and bounded (default at most 64 MiB); no
    transcript sequences are cached.
    The reader is intended for one thread; use a separate reader per worker.
    """

    def __init__(self, database: str | Path, cached_chunks: int = 64) -> None:
        """Open an existing preprocessed database and configure the chunk cache."""
        if cached_chunks < 0:
            raise ValueError("cached_chunks must be nonnegative")
        path = Path(database).expanduser().resolve()
        self.db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
        self.cached_chunks = cached_chunks
        self.chunk_cache: ChunkCache = OrderedDict()
        try:
            metadata = dict(self.db.execute("SELECT key, value FROM build_metadata"))
            if metadata.get("preprocessing_version") != "1":
                raise ValueError("Database is not preprocessed; run --preprocess-only DB")
            if metadata.get("sequence_chunk_size") != str(CHUNK_SIZE):
                raise ValueError("Unsupported sequence chunk size")
            if metadata.get("transcript_payload_codec") not in {None, TRANSCRIPT_PAYLOAD_CODEC}:
                raise ValueError("Unsupported transcript payload codec")
            self.metadata = metadata
        except BaseException:
            self.db.close()
            raise

    def transcript_error_summary(self) -> dict[str, TranscriptErrorSummary]:
        """Read grouped failures with example IDs without rebuilding the reference.

        Prepared summaries are cheap metadata reads. If one is absent, inspect
        only stored error payloads; usable transcript payloads are never loaded.
        """
        if "preprocessing_errors" in self.metadata:
            return cast(
                dict[str, TranscriptErrorSummary], json.loads(self.metadata["preprocessing_errors"])
            )
        errors: dict[str, TranscriptErrorSummary] = {}
        for transcript_id, encoded in self.db.execute(
            "SELECT transcript_id, payload FROM ff_transcripts WHERE status='error' ORDER BY transcript_id"
        ):
            payload = cast("EnsemblError", decode_transcript_payload(encoded))
            message = payload["error"].removeprefix(transcript_id + ": ")
            _record_transcript_error(errors, transcript_id, message)
        return dict(sorted(errors.items(), key=lambda item: (-item[1]["count"], item[0])))

    def get_transcript(
        self, transcript_id: str, *, include_sequence: bool = True
    ) -> ProteinFeatureResponse:
        """Unknown IDs, version mismatches and unsupported models return errors."""
        match = re.fullmatch(r"(ENST\d+)(?:\.(\d+))?", transcript_id)
        if not match:
            return {"error": f"Invalid human Ensembl transcript ID: {transcript_id}"}
        stable_id, requested_version = match.groups()
        row = self.db.execute(
            "SELECT version, status, payload FROM ff_transcripts WHERE transcript_id=?",
            (stable_id,),
        ).fetchone()
        if row is None:
            return {"error": f"Transcript {transcript_id} is absent from this reference release"}
        version, status, encoded = row
        if requested_version is not None and int(requested_version) != version:
            return {
                "error": f"Transcript version mismatch: requested {transcript_id}; stored {stable_id}.{version}"
            }
        if status == "noncoding":
            return {"error": "Transcript is non-coding or has no translation object"}
        payload = decode_transcript_payload(encoded)
        if status == "error":
            return payload
        payload = cast("TranscriptProteinFeatureResult", payload)
        if include_sequence:
            payload["premrna_sequence"] = self.sequence(
                "dna",
                cast(str, payload["chromosome"]),
                payload["transcript_genomic_start"],
                payload["transcript_genomic_end"],
                payload["strand"],
            )
        return payload

    def sequence(
        self, kind: str, sequence_id: str, start: int = 1, end: int | None = None, strand: int = 1
    ) -> str:
        """Read a sequence interval through this reader's shared bounded cache."""
        return fetch_sequence(
            self.db,
            kind,
            sequence_id,
            start,
            end,
            strand,
            _chunk_cache=self.chunk_cache,
            _cache_limit=self.cached_chunks,
        )

    def get_interpro_annotation(self, accession: str) -> ProteinFeatureAnnotation | None:
        """Return stored metadata for an accession, or None when it is absent."""
        row = self.db.execute(
            "SELECT name, entry_type FROM ff_interpro WHERE interpro_id=?", (accession,)
        ).fetchone()
        if row is None:
            return None
        return {"name": row[0], "entry_type": row[1], "interpro_id": accession}

    def close(self) -> None:
        """Release both the database handle and decompressed genome chunks."""
        self.chunk_cache.clear()
        self.db.close()

    def __enter__(self) -> Self:
        """Use this reader as a context manager."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Close the reader even when annotation raises an exception."""
        self.close()


# Full-build orchestration and the prepare-data command.


@contextmanager
def build_lock(path: Path) -> Iterator[None]:
    """Keep concurrent builders from modifying the same resumable database.

    Advisory locks are released by the OS even after a killed process. Leave the
    small lock file in place so another process cannot lock a different inode.
    """
    with path.open("a+b") as stream:
        try:
            if os.name == "nt":
                import msvcrt

                stream.write(b"\0")
                stream.flush()
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)  # type: ignore[attr-defined]  # Windows-only API
            else:
                import fcntl

                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ValueError(f"Another build is using this output: {path}") from exc
        try:
            yield
        finally:
            if os.name == "nt":
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)  # type: ignore[attr-defined]  # Windows-only API
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


class BuildCheckpoint:
    """Record completed stages only after their SQLite work has committed."""

    def __init__(self, db: sqlite3.Connection, configuration: Mapping[str, object]) -> None:
        self.db = db
        db.execute(
            "CREATE TABLE IF NOT EXISTS build_checkpoints ("
            "stage TEXT PRIMARY KEY, source TEXT NOT NULL, fingerprint TEXT NOT NULL, "
            "count INTEGER NOT NULL, complete INTEGER NOT NULL)"
        )
        fingerprint = json.dumps(configuration, sort_keys=True)
        previous = self.read("configuration")
        if previous and previous["fingerprint"] != fingerprint:
            raise ValueError(
                "Build checkpoint belongs to different source/settings; "
                "use a different --output or remove its .building database"
            )
        if previous is None:
            self.save("configuration", fingerprint=fingerprint)

    def read(self, stage: str) -> CheckpointRecord | None:
        """Read one stage's source provenance and completion marker."""
        row = self.db.execute(
            "SELECT source, fingerprint, count, complete FROM build_checkpoints WHERE stage=?",
            (stage,),
        ).fetchone()
        if row is None:
            return None
        return {
            "source": json.loads(row[0]),
            "fingerprint": row[1],
            "count": row[2],
            "complete": bool(row[3]),
        }

    def save(
        self,
        stage: str,
        *,
        source: Mapping[str, Any] | None = None,
        fingerprint: str = "",
        count: int = 0,
        complete: bool = True,
    ) -> None:
        """Commit a marker with the work it describes; incomplete stages can resume."""
        self.db.execute(
            "INSERT OR REPLACE INTO build_checkpoints VALUES (?, ?, ?, ?, ?)",
            (stage, json.dumps(source or {}, sort_keys=True), fingerprint, count, int(complete)),
        )
        self.db.commit()


def build(args: argparse.Namespace) -> Path:
    """Resume committed imports in a staging database, then publish atomically.

    Failures preserve the staging database and downloads. Repeating the same
    command skips completed stages; no incomplete reference replaces the output.
    """
    from .uniprot import UNIPROT_HUMAN_URL

    cache = args.cache_dir.expanduser().resolve()
    base = args.base_url.rstrip("/") + "/"
    # Resolution, metadata/schema, tables, FASTA, preprocessing, analysis, check, publication.
    steps = BuildSteps(1 + 3 + len(TABLES) + len(SEQUENCE_TYPES) + 4)
    with steps.step("Resolve Ensembl release and human GRCh38 core directory"):
        release, root, core = resolve_core(base, args.release)
    LOG.info("Resolved Ensembl release %s; core database %s", release, core)
    output = (
        (args.output or cache / args.species / ASSEMBLY / f"release-{release}" / "ensembl.sqlite")
        .expanduser()
        .resolve()
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.force:
        raise FileExistsError(f"Output exists: {output}. Use --force to replace it.")
    downloads = cache / "downloads" / hashlib.sha256(base.encode()).hexdigest()[:12] / core
    interpro_cache = cache / "metadata" / "interpro"
    LOG.info("Final reference database: %s", output)
    for group in ("core", *SEQUENCE_TYPES):
        LOG.info("%s download cache: %s", group, downloads / group)
    LOG.info("InterPro download cache: %s", interpro_cache)
    LOG.info("Historical InterPro entry-list cache: %s", interpro_cache / "releases")
    LOG.info(
        "Partial downloads: <download-cache>/<filename>.<random>.part; checksums: <filename>.sha256"
    )
    temporary = output.with_name(output.name + ".building")
    lock_path = output.with_name(output.name + ".prepare.lock")
    LOG.info("Resumable build checkpoint: %s; SQLite journal: %s-journal", temporary, temporary)
    LOG.info("Build lock: %s", lock_path)
    source_rows = []
    fetch_interpro = args.interpro_entries is None

    def get(url: str, group: str = "core", *, directory: Path | None = None) -> Path:
        """Download/cache one source and retain its checksum/size for provenance."""
        path, digest = download(url, directory if directory is not None else downloads / group)
        source_rows.append((url, digest, path.stat().st_size))
        return path

    core_url = root + f"mysql/{core}/"
    LOG.info("Temporary reference database: %s; rerun the same command to resume", temporary)
    try:
        with steps.step("Prepare InterPro metadata"):
            if fetch_interpro:
                interpro_entries = get(INTERPRO_ENTRIES_URL, directory=interpro_cache)
            else:
                interpro_entries = args.interpro_entries
                LOG.info("Using local InterPro entries: %s", interpro_entries)
        with steps.step("Prepare reviewed human UniProt features"):
            if args.uniprot_features is None:
                uniprot_features = get(UNIPROT_HUMAN_URL, directory=cache / "metadata" / "uniprot")
                uniprot_source = {"url": UNIPROT_HUMAN_URL, "sha256": source_rows[-1][1]}
            else:
                uniprot_features = args.uniprot_features
                uniprot_source = {
                    "path": str(uniprot_features),
                    "sha256": sha256_file(uniprot_features),
                }
                LOG.info("Using local UniProt features: %s", uniprot_features)
        with steps.step("Download and read Ensembl table definitions"):
            schema = parse_schema(get(core_url + core + ".sql.gz"))
            absent = set(TABLES) - set(schema)
            if absent:
                raise ValueError(f"Required tables absent from source schema: {sorted(absent)}")
        with build_lock(lock_path), closing(sqlite3.connect(temporary)) as db:
            if output.exists() and not args.force:
                raise FileExistsError(f"Output appeared during the build: {output}")
            # Apply build-only settings to the unpublished database. The 64 MiB
            # page cache reduces repeated reads without caching the entire input.
            db.execute("PRAGMA journal_mode=DELETE")
            db.execute("PRAGMA synchronous=NORMAL")
            db.execute("PRAGMA cache_size=-65536")
            db.execute("PRAGMA user_version=1")
            checkpoint = BuildCheckpoint(
                db,
                {
                    "checkpoint_version": 1,
                    "base_url": base,
                    "release": release,
                    "species": args.species,
                    "core_database": core,
                    "assembly": ASSEMBLY,
                    "schema_sha256": source_rows[-1][1],
                    "sequence_chunk_size": CHUNK_SIZE,
                    "tables": TABLES,
                    "sequence_types": SEQUENCE_TYPES,
                },
            )

            def reuse(stage: str) -> int | None:
                """Reuse a completed import and its recorded source provenance."""
                saved = checkpoint.read(stage)
                if saved and saved["complete"]:
                    source = saved["source"]
                    source_rows.append((source["url"], source["sha256"], source["bytes"]))
                    LOG.info("Reusing checkpoint for %s: %s records", stage, f"{saved['count']:,}")
                    return saved["count"]
                return None

            def last_source() -> SourceRecord:
                """Describe the most recently verified download."""
                url, digest, size = source_rows[-1]
                return {"url": url, "sha256": digest, "bytes": size}

            counts = {}
            panther_digest = None
            for table in TABLES:
                with steps.step(f"Download, import and index {table}"):
                    count = reuse("table:" + table)
                    if count is None:
                        # A failed import may have created its table but has no
                        # completion marker. Recreate only that unfinished table.
                        db.execute(f"DROP TABLE IF EXISTS {sql_name('ensembl_' + table)}")
                        count = import_table(
                            db, table, schema[table], get(core_url + table + ".txt.gz")
                        )
                        checkpoint.save("table:" + table, source=last_source(), count=count)
                    counts[table] = count
                    if table == "analysis" and args.panther_classifications is None:
                        version = panther_release(protein_analyses(db))
                        if version:
                            LOG.info(
                                "Check/download PANTHER %s metadata before importing genome sequences",
                                version,
                            )
                            _, panther_digest, _ = download_panther_classifications(
                                version, cache / "metadata" / "panther"
                            )
            # Guard against accidentally mixing source database releases.
            versions = {
                str(row[0])
                for row in db.execute(
                    "SELECT meta_value FROM ensembl_meta WHERE meta_key='schema_version'"
                )
            }
            if versions and versions != {str(release)}:
                raise ValueError(f"Core schema_version {versions} does not match release {release}")
            if not checkpoint.read("sequence_tables"):
                # SQLite executes CREATE TABLE outside implicit DML transactions;
                # recover a process killed partway through schema creation.
                db.execute("DROP TABLE IF EXISTS sequence_chunks")
                db.execute("DROP TABLE IF EXISTS sequences")
                create_sequence_tables(db)
                checkpoint.save("sequence_tables")
            for kind in SEQUENCE_TYPES:
                with steps.step(f"Download and import {kind} FASTA sequences"):
                    count = reuse("fasta:" + kind)
                    if count is not None:
                        counts[kind + "_sequences"] = count
                        continue
                    directory = root + f"fasta/{args.species}/{kind}/"
                    suffix = ".dna.toplevel.fa.gz" if kind == "dna" else f".{kind}.all.fa.gz"
                    candidates = sorted(
                        name for name in listing(directory) if name.endswith(suffix)
                    )
                    if len(candidates) != 1:
                        raise ValueError(
                            f"Expected one {suffix} file at {directory}; found {candidates}"
                        )
                    path = get(directory + candidates[0], kind)
                    source = last_source()
                    saved = checkpoint.read("fasta:" + kind)
                    if saved is None or saved["source"] != source:
                        # Never reuse records after the source bytes change.
                        db.execute("DELETE FROM sequence_chunks WHERE kind=?", (kind,))
                        db.execute("DELETE FROM sequences WHERE kind=?", (kind,))
                    checkpoint.save("fasta:" + kind, source=source, complete=False)
                    count = import_fasta(db, kind, path, resume=True)
                    counts[kind + "_sequences"] = count
                    if kind == "dna":
                        validate_dna_regions(db)
                    checkpoint.save("fasta:" + kind, source=source, count=count)
            db.execute(
                "CREATE TABLE IF NOT EXISTS build_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            metadata = {
                "format_version": "1",
                "release": str(release),
                "species": args.species,
                "core_database": core,
                "assembly": ASSEMBLY,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "sequence_chunk_size": str(CHUNK_SIZE),
                "sequence_codec": "zlib",
                "tables": json.dumps(TABLES),
                "sequence_types": json.dumps(SEQUENCE_TYPES),
                "counts": json.dumps(counts),
            }
            db.executemany("INSERT OR REPLACE INTO build_metadata VALUES (?, ?)", metadata.items())
            db.execute(
                "CREATE TABLE IF NOT EXISTS source_files (url TEXT PRIMARY KEY, sha256 TEXT, bytes INTEGER)"
            )
            db.executemany("INSERT OR REPLACE INTO source_files VALUES (?, ?, ?)", source_rows)
            db.commit()
            with steps.step("Preprocess transcripts, domains and splice sites"):
                # If preprocessing completed before a later failure, preserve it
                # unless its metadata input or implementation has changed.
                fingerprint = json.dumps(
                    {
                        "implementation": sha256_file(Path(__file__)),
                        "uniprot_implementation": sha256_file(
                            Path(__file__).with_name("uniprot.py")
                        ),
                        "ensembl_implementation": sha256_file(
                            Path(__file__).with_name("ensembl.py")
                        ),
                        "uniprot_sha256": sha256_file(uniprot_features),
                        "interpro_sha256": sha256_file(interpro_entries),
                        "fetch_interpro": fetch_interpro,
                        "panther_sha256": (
                            sha256_file(args.panther_classifications)
                            if args.panther_classifications
                            else panther_digest
                        ),
                    },
                    sort_keys=True,
                )
                saved = checkpoint.read("preprocessing")
                if saved and saved["complete"] and saved["fingerprint"] == fingerprint:
                    LOG.info("Reusing completed transcript/domain preprocessing")
                else:
                    db.execute("DELETE FROM build_checkpoints WHERE stage='analyze'")
                    checkpoint.save("preprocessing", fingerprint=fingerprint, complete=False)
                    preprocess_reference(
                        db,
                        interpro_entries,
                        interpro_archive_cache=interpro_cache if fetch_interpro else None,
                        panther_classifications=args.panther_classifications,
                        panther_cache=cache / "metadata" / "panther",
                        uniprot_features=uniprot_features,
                        uniprot_source=uniprot_source,
                    )
                    checkpoint.save("preprocessing", fingerprint=fingerprint)
            with steps.step("Analyze database indexes"):
                if checkpoint.read("analyze"):
                    LOG.info("Reusing completed index analysis")
                else:
                    sqlite_phase(db, "ANALYZE", "Analyze database")
                    checkpoint.save("analyze")
            with steps.step("Validate SQLite database integrity"):
                result = sqlite_phase(db, "PRAGMA integrity_check", "Check database")
                if result != [("ok",)]:
                    raise ValueError(f"SQLite integrity check failed: {result[:5]}")
            db.commit()
            # Close SQLite before publication, while still holding the build
            # lock. This also permits atomic replacement on Windows.
            db.close()
            with steps.step("Publish completed reference database"):
                if output.exists() and not args.force:
                    raise FileExistsError(f"Output appeared during the build: {output}")
                os.replace(temporary, output)
    except BaseException:
        if temporary.exists():
            LOG.error(
                "Build checkpoint preserved at %s; rerun the same command to resume", temporary
            )
        raise
    LOG.info("Created %s (%.1f MiB)", output, output.stat().st_size / CHUNK_SIZE)
    steps.complete()
    return output


def main(argv: Sequence[str] | None = None) -> int:
    """Validate CLI options, then run a full build or derived-table regeneration."""
    parser = argparse.ArgumentParser(
        prog="fusion-function prepare-data",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--release", type=int, help="Ensembl release; default: latest numbered FTP release"
    )
    parser.add_argument(
        "--species",
        choices=[SPECIES],
        help="Species component of the output path; this script supports human data",
    )
    parser.add_argument(
        "--output", type=Path, help="SQLite file; default: release-specific file in cache"
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=default_cache_dir(),
        help="Download/cache directory; defaults to FUSION_FUNCTION_CACHEDIR or OS user cache",
    )
    parser.add_argument("--base-url", help="Advanced: FTP HTTP(S) mirror root, ending at /pub/")
    parser.add_argument(
        "--from-source",
        action="store_true",
        default=None,
        help="Build from Ensembl FTP instead of downloading a compatible prebuilt reference",
    )
    parser.add_argument(
        "--reference-catalog",
        help="Prebuilt catalog: local JSON or HTTPS URL; default: maintained catalog",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=None,
        help="Replace output after successful installation or build",
    )
    parser.add_argument(
        "--preprocess-only",
        type=Path,
        metavar="DB",
        help="Regenerate derived tables in a full source-built database",
    )
    parser.add_argument(
        "--interpro-entries",
        type=Path,
        metavar="TSV",
        help="Optional local InterPro entry.list TSV for canonical entry names/types",
    )
    parser.add_argument(
        "--panther-classifications",
        type=Path,
        metavar="TSV",
        help="Optional local PANTHER human classification TSV; default: Ensembl's PANTHER release",
    )
    parser.add_argument(
        "--uniprot-features",
        type=Path,
        metavar="XML",
        help="Local UniProt XML (optionally gzip); default: reviewed human bulk download",
    )
    args = parser.parse_args(argv)
    if args.release is not None and args.release <= 0:
        parser.error("--release must be positive")
    if args.preprocess_only:
        # Reprocessing uses the database's source release/assembly, so new-build
        # options cannot change where its existing coordinates came from.
        incompatible = [
            "--" + name.replace("_", "-")
            for name in (
                "release",
                "species",
                "output",
                "base_url",
                "force",
                "from_source",
                "reference_catalog",
            )
            if getattr(args, name) is not None
        ]
        if incompatible:
            parser.error("--preprocess-only cannot be combined with " + ", ".join(incompatible))
    source_options = args.from_source or any(
        getattr(args, name) is not None
        for name in ("base_url", "interpro_entries", "panther_classifications", "uniprot_features")
    )
    if args.reference_catalog and source_options:
        parser.error("--reference-catalog cannot be combined with source-build options")
    if args.interpro_entries:
        args.interpro_entries = args.interpro_entries.expanduser().resolve()
        if not args.interpro_entries.is_file():
            parser.error(f"InterPro entry file does not exist: {args.interpro_entries}")
    if args.panther_classifications:
        args.panther_classifications = args.panther_classifications.expanduser().resolve()
        if not args.panther_classifications.is_file():
            parser.error(
                f"PANTHER classification file does not exist: {args.panther_classifications}"
            )
    if args.uniprot_features:
        args.uniprot_features = args.uniprot_features.expanduser().resolve()
        if not args.uniprot_features.is_file():
            parser.error(f"UniProt feature file does not exist: {args.uniprot_features}")
    args.species = args.species or SPECIES
    args.base_url = args.base_url or BASE_URL
    args.force = bool(args.force)
    parsed_base = urllib.parse.urlparse(args.base_url)
    if (
        parsed_base.scheme not in {"http", "https"}
        or not parsed_base.netloc
        or parsed_base.query
        or parsed_base.fragment
    ):
        parser.error("--base-url must be an HTTP(S) mirror URL without a query or fragment")
    args.cache_dir = args.cache_dir.expanduser().resolve()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    LOG.info("Cache root: %s", args.cache_dir)
    LOG.info("Ensembl download root: %s", args.cache_dir / "downloads")
    LOG.info("InterPro cache: %s", args.cache_dir / "metadata" / "interpro")
    LOG.info(
        "Historical InterPro entry-list cache: %s",
        args.cache_dir / "metadata" / "interpro" / "releases",
    )
    LOG.info("UniProt feature cache: %s", args.cache_dir / "metadata" / "uniprot")
    LOG.info("PANTHER human classification cache: %s", args.cache_dir / "metadata" / "panther")
    LOG.info("Prebuilt download cache: %s", args.cache_dir / "prebuilt")
    LOG.info("Source partial downloads: <filename>.<random>.part beside cached files")
    LOG.info("Prebuilt partial downloads: <cache>/prebuilt/<sha256>/ensembl.sqlite.gz.part")
    if args.interpro_entries:
        LOG.info("Local InterPro entries: %s", args.interpro_entries)
    if args.panther_classifications:
        LOG.info("Local PANTHER classifications: %s", args.panther_classifications)
    with logging_redirect_tqdm():
        try:
            if args.preprocess_only:
                output = args.preprocess_only.expanduser().resolve()
                LOG.info("Reprocessing existing reference in place: %s", output)
                LOG.info(
                    "SQLite transaction files, if used: %s-journal; %s-wal; %s-shm",
                    output,
                    output,
                    output,
                )
                steps = BuildSteps(4)
                with (
                    build_lock(output.with_name(output.name + ".prepare.lock")),
                    closing(sqlite3.connect(output.as_uri() + "?mode=rw", uri=True)) as db,
                ):
                    with steps.step("Validate existing reference database"):
                        validate_reference_source(db)
                    with steps.step("Prepare InterPro metadata"):
                        if args.interpro_entries is None:
                            archive_cache = args.cache_dir / "metadata" / "interpro"
                            interpro_entries, _ = download(INTERPRO_ENTRIES_URL, archive_cache)
                        else:
                            interpro_entries = args.interpro_entries
                            archive_cache = None
                            LOG.info("Using local InterPro entries: %s", interpro_entries)
                    with steps.step("Prepare reviewed human UniProt features"):
                        from .uniprot import UNIPROT_HUMAN_URL

                        if args.uniprot_features is None:
                            uniprot_features, digest = download(
                                UNIPROT_HUMAN_URL, args.cache_dir / "metadata" / "uniprot"
                            )
                            uniprot_source = {"url": UNIPROT_HUMAN_URL, "sha256": digest}
                        else:
                            uniprot_features = args.uniprot_features
                            uniprot_source = {
                                "path": str(uniprot_features),
                                "sha256": sha256_file(uniprot_features),
                            }
                            LOG.info("Using local UniProt features: %s", uniprot_features)
                    with steps.step("Preprocess transcripts, domains and splice sites"):
                        preprocess_reference(
                            db,
                            interpro_entries,
                            interpro_archive_cache=archive_cache,
                            panther_classifications=args.panther_classifications,
                            panther_cache=args.cache_dir / "metadata" / "panther",
                            uniprot_features=uniprot_features,
                            uniprot_source=uniprot_source,
                        )
                steps.complete()
            else:
                LOG.info(
                    "Reference destination: %s",
                    args.output.expanduser().resolve()
                    if args.output
                    else args.cache_dir
                    / args.species
                    / ASSEMBLY
                    / f"release-{args.release or '<latest>'}"
                    / "ensembl.sqlite",
                )
                output = None
                if not source_options:
                    from .prebuilt import install_reference

                    output = install_reference(
                        release=args.release,
                        cache_dir=args.cache_dir,
                        output=args.output,
                        force=args.force,
                        catalog=args.reference_catalog,
                    )
                else:
                    LOG.info("Source-build options selected; skipping the prebuilt catalog")
                if output is None:
                    output = build(args)
        except (
            OSError,
            ValueError,
            sqlite3.Error,
            urllib.error.URLError,
            EOFError,
            ParseError,
        ) as exc:
            LOG.error("Build failed: %s", exc)
            return 1
        except KeyboardInterrupt:
            LOG.error("Interrupted; no incomplete database was published")
            return 130
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
