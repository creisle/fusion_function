#!/usr/bin/env python3
"""Benchmark an existing fusion-function reference without altering it.

Examples (run in the package's Python environment):
    python scripts/benchmark_fusion_function.py --release 116 --output benchmark-116.json
    python scripts/benchmark_fusion_function.py --database /path/ensembl.sqlite
    python scripts/benchmark_fusion_function.py --cases sv_cases.json --sample 100 --repeat 3

--cases accepts a JSON list of annotate_fusion_domains inputs, with optional
"name" fields. Reference/database/release arguments belong on this script's CLI.
Default cases use the package's GRCh38 regression breakpoints: BCR::ABL1,
EML4::ALK, KIF5B::RET, EWSR1::FLI1 and SLC45A2::AMACR.

First pass means first access in THIS PROCESS, not a guaranteed cold disk cache.
The script never drops OS caches. Repeat passes reuse the reader and its bounded
sequence cache, just like a batch annotation job. Diverse transcript profiles
run AFTER fusion timings so profiling does not pre-warm fusion inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import sqlite3
import statistics
import sys
import time
from collections import defaultdict
from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, TypeVar, cast

from tqdm import tqdm

from fusion_function import ReferenceDatabase, annotate_fusion_domains
from fusion_function.data import decode_transcript_payload

T = TypeVar("T")
Record = dict[str, Any]


def default_cases() -> list[Record]:
    """Use fixed regression inputs, avoiding transcript reads during setup."""
    presets = [
        ("BCR::ABL1", "ENST00000305877", "ENST00000318560", "22:23290413", "9:130854064"),
        ("EML4::ALK", "ENST00000318522", "ENST00000389048", "2:42295516", "2:29223528"),
        ("KIF5B::RET", "ENST00000302418", "ENST00000355710", "10:32028428", "10:43116584"),
        ("EWSR1::FLI1", "ENST00000397938", "ENST00000527786", "22:29287870", "11:128793694"),
        ("SLC45A2::AMACR", "ENST00000296589", "ENST00000335606", "5:33973126", "5:34006836"),
    ]
    return [
        dict(
            name=name,
            transcript1_id=t1,
            transcript2_id=t2,
            breakpoint1=b1,
            breakpoint2=b2,
            gene1_terminus="N",
            gene2_terminus="C",
        )
        for name, t1, t2, b1, b2 in presets
    ]


def load_cases(path: Path | None) -> list[Record]:
    """Reject misspelled/unsupported inputs instead of silently ignoring them."""
    cases = default_cases() if path is None else json.loads(path.read_text())
    allowed = {
        "name",
        "transcript1_id",
        "transcript2_id",
        "breakpoint1",
        "breakpoint2",
        "gene1_terminus",
        "gene2_terminus",
        "inserted_sequence",
    }
    if not isinstance(cases, list) or not cases:
        raise ValueError("Cases must be a nonempty JSON list")
    for number, case in enumerate(cases, 1):
        if not isinstance(case, dict) or set(case) - allowed:
            raise ValueError(f"Case {number}: expected an object using only {sorted(allowed)}")
        if any(value is not None and not isinstance(value, str) for value in case.values()):
            raise ValueError(f"Case {number}: values must be strings or null")
        if not case.get("transcript1_id") and not case.get("transcript2_id"):
            raise ValueError(f"Case {number}: at least one transcript is required")
        for side in (1, 2):
            if case.get(f"transcript{side}_id") and (
                not case.get(f"breakpoint{side}")
                or case.get(f"gene{side}_terminus") not in {"N", "C"}
            ):
                raise ValueError(f"Case {number}: known partner {side} needs a breakpoint and N/C")
        case["name"] = case.get("name") or f"case-{number}"
    return cases


def sample_transcripts(reference: ReferenceDatabase, count: int, seed: int) -> list[str]:
    """Uniform reservoir sample over ready IDs; do not read their JSON payloads."""
    if count == 0:
        return []
    rng = random.Random(seed)
    selected: list[str] = []
    query = "SELECT transcript_id FROM ff_transcripts WHERE status='ready' ORDER BY transcript_id"
    for number, (transcript_id,) in enumerate(reference.db.execute(query), 1):
        if number <= count:
            selected.append(transcript_id)
        else:
            index = rng.randrange(number)
            if index < count:
                selected[index] = transcript_id
    rng.shuffle(selected)
    return selected


def measure(
    rows: list[Record], operation: str, phase: str, name: str, call: Callable[[], T]
) -> tuple[T | None, Record]:
    """Time only the operation, excluding progress updates and result hashing."""
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    result: T | None = None
    error: str | None = None
    try:
        result = call()
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    wall = time.perf_counter() - wall_start
    cpu = time.process_time() - cpu_start
    row = dict(
        operation=operation, phase=phase, name=name, wall_seconds=wall, cpu_seconds=cpu, error=error
    )
    rows.append(row)
    return result, row


def annotate(reference: ReferenceDatabase, case: Record) -> Record:
    """Call the production API with explicit supported arguments."""
    return dict(
        annotate_fusion_domains(
            transcript1_id=case.get("transcript1_id"),
            transcript2_id=case.get("transcript2_id"),
            breakpoint1=case.get("breakpoint1"),
            breakpoint2=case.get("breakpoint2"),
            gene1_terminus=case.get("gene1_terminus"),
            gene2_terminus=case.get("gene2_terminus"),
            inserted_sequence=case.get("inserted_sequence") or "",
            reference=reference,
        )
    )


def benchmark_fusions(
    reference: ReferenceDatabase, cases: list[Record], phases: list[str], rows: list[Record]
) -> None:
    """Measure full calls, flag returned errors and inconsistent repeat results."""
    fingerprints: dict[int, str] = {}
    for phase in phases:
        with tqdm(cases, desc=f"[3/5] Fusion {phase}", unit="SV") as bar:
            for number, case in enumerate(bar):
                bar.set_postfix_str(case["name"], refresh=False)
                result, row = measure(
                    rows, "fusion", phase, case["name"], lambda: annotate(reference, case)
                )
                row["case_index"] = number
                if result is not None:
                    row["error"] = result.get("error")
                    if not row["error"]:
                        fingerprint = hashlib.sha256(
                            json.dumps(result, sort_keys=True).encode()
                        ).hexdigest()
                        previous = fingerprints.setdefault(number, fingerprint)
                        row.update(
                            result_sha256=fingerprint,
                            frame_status=result.get("frame_status"),
                            features=len(result.get("domains", [])),
                        )
                        if previous != fingerprint:
                            row["error"] = "Annotation changed between passes"
                if row["error"]:
                    tqdm.write(f"  {case['name']}: {row['error']}", file=sys.stderr)


def benchmark_transcripts(
    reference: ReferenceDatabase, ids: list[str], phases: list[str], rows: list[Record]
) -> None:
    """Separate SQL retrieval, model decoding and production sequence retrieval.

    Model decoding includes decompression for compressed transcript records.

    Sequence retrieval includes SQL, decompression, slicing and strand handling.
    Full pre-mRNA spans introns, so it can require multiple 1 MiB genome chunks.
    These are component timings, not additional full annotation measurements.
    """
    for phase in phases:
        with tqdm(ids, desc=f"[4/5] Transcripts {phase}", unit="transcript") as bar:
            for transcript_id in bar:
                # Let tqdm throttle rendering. Forced redraws can dominate a
                # short repeat pass over a remote terminal and are not API work.
                bar.set_postfix_str(f"{transcript_id}: SQL", refresh=False)
                stored, sql_row = measure(
                    rows,
                    "transcript_sql",
                    phase,
                    transcript_id,
                    lambda: reference.db.execute(
                        "SELECT payload FROM ff_transcripts WHERE transcript_id=?", (transcript_id,)
                    ).fetchone(),
                )
                if stored is None:
                    sql_row["error"] = sql_row["error"] or "Transcript not found"
                    continue
                sql_row["payload_bytes"] = len(
                    stored[0].encode() if isinstance(stored[0], str) else stored[0]
                )
                bar.set_postfix_str(f"{transcript_id}: JSON", refresh=False)
                payload, parse_row = measure(
                    rows,
                    "transcript_json",
                    phase,
                    transcript_id,
                    lambda: cast(Record, decode_transcript_payload(stored[0])),
                )
                if payload is None:
                    continue
                if "error" in payload:
                    parse_row["error"] = payload["error"]
                    continue
                bar.set_postfix_str(f"{transcript_id}: sequence", refresh=False)
                sequence, sequence_row = measure(
                    rows,
                    "transcript_sequence",
                    phase,
                    transcript_id,
                    lambda: reference.sequence(
                        "dna",
                        cast(str, payload["chromosome"]),
                        payload["transcript_genomic_start"],
                        payload["transcript_genomic_end"],
                        payload["strand"],
                    ),
                )
                sequence_row.update(
                    premrna_bases=payload["premrna_length"],
                    features=len(payload["protein_features"]),
                )
                if sequence is not None and len(sequence) != payload["premrna_length"]:
                    sequence_row["error"] = "Sequence length differs from transcript model"


def summarize(rows: list[Record]) -> list[Record]:
    """Report successful latency distributions separately from failed calls."""
    groups: dict[tuple[str, str], list[Record]] = defaultdict(list)
    for row in rows:
        groups[row["operation"], row["phase"]].append(row)
    summaries = []
    for (operation, phase), values in groups.items():
        wall = sorted(row["wall_seconds"] for row in values if not row["error"])
        summaries.append(
            dict(
                operation=operation,
                phase=phase,
                successful=len(wall),
                errors=len(values) - len(wall),
                wall_total_seconds=sum(row["wall_seconds"] for row in values),
                cpu_total_seconds=sum(row["cpu_seconds"] for row in values),
                median_seconds=statistics.median(wall) if wall else None,
                p95_seconds=wall[max(0, (95 * len(wall) + 99) // 100 - 1)] if wall else None,
                max_seconds=max(wall) if wall else None,
            )
        )
    return summaries


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--database", type=Path, help="Existing DB; default: package cache discovery"
    )
    parser.add_argument("--release", type=int, help="Select an installed Ensembl release")
    parser.add_argument("--cases", type=Path, help="JSON list of actual SV API inputs")
    parser.add_argument(
        "--sample", type=int, default=100, help="Random ready transcript profiles (0 disables)"
    )
    parser.add_argument("--repeat", type=int, default=3, help="Repeat passes after the first pass")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--cached-chunks", type=int, default=64, help="Production default: 64 chunks"
    )
    parser.add_argument(
        "--output", type=Path, help="New JSON result file; existing files are not overwritten"
    )
    args = parser.parse_args(argv)
    if min(args.sample, args.repeat, args.cached_chunks) < 0:
        parser.error("sample, repeat and cached-chunks must be nonnegative")
    try:
        cases = load_cases(args.cases)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    output = (
        (
            args.output
            or Path(datetime.now(UTC).strftime("benchmark-fusion-function-%Y%m%dT%H%M%S%fZ.json"))
        )
        .expanduser()
        .resolve()
    )
    if args.database and output == args.database.expanduser().resolve():
        parser.error("Output must not be the reference database")
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        stream = output.open("x")
    except OSError as exc:
        parser.error(str(exc))
    rows: list[Record] = []
    try:
        package_version = version("fusion-function")
    except PackageNotFoundError:
        package_version = "source checkout"
    report: Record = dict(
        status="running",
        started_utc=datetime.now(UTC).isoformat(),
        python=sys.version,
        platform=platform.platform(),
        sqlite_version=sqlite3.sqlite_version,
        package_version=package_version,
        seed=args.seed,
        cached_chunks=args.cached_chunks,
        repeats=args.repeat,
        cases=cases,
        measurements=rows,
        cache_note="First pass is not guaranteed OS-cold; repeat passes reuse bounded caches.",
    )
    exit_code = 0
    started = time.perf_counter()
    try:
        print(f"Results: {output}", flush=True)
        print("[1/5] Open existing reference (read-only)", flush=True)
        reference, opened = measure(
            rows,
            "open_reference",
            "first",
            "reference",
            lambda: ReferenceDatabase(
                args.database, release=args.release, cached_chunks=args.cached_chunks
            ),
        )
        if reference is None:
            raise RuntimeError(opened["error"])
        with reference:
            report.update(
                database=str(reference.path),
                database_bytes=reference.path.stat().st_size,
                release=reference.metadata.get("release"),
                preprocessed_utc=reference.metadata.get("preprocessed_utc"),
            )
            report["sqlite_settings"] = {
                setting: reference.db.execute(f"PRAGMA {setting}").fetchone()[0]
                for setting in ("page_size", "cache_size", "journal_mode", "mmap_size")
            }
            transcript_schema = reference.db.execute(
                "SELECT sql FROM sqlite_master WHERE type='table' AND name='ff_transcripts'"
            ).fetchone()[0]
            report["transcript_storage"] = {
                "payload_codec": reference.metadata.get("transcript_payload_codec", "json"),
                "table_layout": (
                    "without_rowid" if "WITHOUT ROWID" in transcript_schema.upper() else "rowid"
                ),
            }
            print(f"Database: {reference.path} (release {report['release']})", flush=True)
            print(
                "Transcript storage: "
                f"{report['transcript_storage']['payload_codec']}, "
                f"{report['transcript_storage']['table_layout']}",
                flush=True,
            )
            print(
                "First pass is NOT guaranteed cold: OS/storage caches are left intact.", flush=True
            )
            print("[2/5] Select random ready IDs using metadata only", flush=True)
            ids, selected = measure(
                rows,
                "select_ids",
                "setup",
                "sample",
                lambda: sample_transcripts(reference, args.sample, args.seed),
            )
            if ids is None:
                raise RuntimeError(selected["error"])
            report["sampled_transcripts"] = ids
            print(
                f"Workload: {len(cases)} SVs; {len(ids)} transcript profiles; "
                f"{args.repeat + 1} passes",
                flush=True,
            )
            phases = ["first"] + [f"repeat-{number}" for number in range(1, args.repeat + 1)]
            benchmark_fusions(reference, cases, phases, rows)
            benchmark_transcripts(reference, ids, phases, rows)
        report["status"] = "completed"
        exit_code = int(any(row["error"] for row in rows))
    except KeyboardInterrupt:
        report["status"], exit_code = "interrupted", 130
        print("\nInterrupted; saving completed measurements.", file=sys.stderr, flush=True)
    except Exception as exc:
        report["status"], report["error"], exit_code = "failed", str(exc), 1
        print(f"Benchmark failed: {exc}", file=sys.stderr, flush=True)
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        report["summary"] = summarize(rows)
        with stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
        print("[5/5] Timing summary (seconds; latency statistics exclude failed calls)", flush=True)
        print(
            f"{'operation':24} {'pass':10} {'ok/err':8} {'wall':>9} {'CPU':>9} "
            f"{'median':>9} {'p95':>9} {'max':>9}",
            flush=True,
        )
        for item in report["summary"]:
            latencies = [
                f"{item[key]:.4f}" if item[key] is not None else "NA"
                for key in ("median_seconds", "p95_seconds", "max_seconds")
            ]
            print(
                f"{item['operation']:24} {item['phase']:10} "
                f"{item['successful']}/{item['errors']:<6} "
                f"{item['wall_total_seconds']:9.4f} {item['cpu_total_seconds']:9.4f} "
                + " ".join(f"{value:>9}" for value in latencies),
                flush=True,
            )
        slowest = sorted(
            (
                row
                for row in rows
                if row["operation"] in {"fusion", "transcript_sql", "transcript_sequence"}
            ),
            key=lambda row: row["wall_seconds"],
            reverse=True,
        )[:5]
        print("Slowest completed operations:", flush=True)
        for row in slowest:
            print(
                f"  {row['wall_seconds']:.4f}s wall / {row['cpu_seconds']:.4f}s CPU: "
                f"{row['operation']} {row['phase']} {row['name']}"
                + (f" ERROR: {row['error']}" if row["error"] else ""),
                flush=True,
            )
        print(f"Saved: {output}", flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
