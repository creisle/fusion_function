"""Reference configuration and read-only access to prepared human databases."""

from __future__ import annotations

import os
import re
import threading
from pathlib import Path

from .data import ASSEMBLY, SPECIES, ReferenceReader, default_cache_dir


def resolve_database(
    database: str | Path | None = None,
    *,
    release: int | None = None,
    cache_dir: str | Path | None = None,
) -> Path:
    """Choose explicit path, FUSION_FUNCTION_DB, or newest prepared local release."""
    if release is not None and (not isinstance(release, int) or release <= 0):
        raise ValueError("release must be a positive Ensembl release number")
    explicit = database or os.environ.get("FUSION_FUNCTION_DB")
    if explicit:
        path = Path(explicit).expanduser().resolve()
    else:
        root = Path(cache_dir or default_cache_dir()).expanduser() / SPECIES / ASSEMBLY
        if release is not None:
            path = root / f"release-{release}" / "ensembl.sqlite"
        else:
            candidates = []
            if root.is_dir():
                for directory in root.iterdir():
                    match = re.fullmatch(r"release-(\d+)", directory.name)
                    if match and (directory / "ensembl.sqlite").is_file():
                        candidates.append((int(match[1]), directory / "ensembl.sqlite"))
            if not candidates:
                raise FileNotFoundError(
                    f"No prepared human reference in {root}. Run 'fusion-function prepare-data' "
                    "once, or set FUSION_FUNCTION_DB to an existing preprocessed database."
                )
            path = max(candidates)[1]
        path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Reference database not found: {path}. Run 'fusion-function prepare-data'."
        )
    return path


class ReferenceDatabase(ReferenceReader):
    """Reusable per-thread reference reader; close explicitly or use a context manager."""

    def __init__(
        self,
        database: str | Path | None = None,
        *,
        release: int | None = None,
        cache_dir: str | Path | None = None,
        cached_chunks: int = 64,
    ) -> None:
        self.path = resolve_database(database, release=release, cache_dir=cache_dir)
        super().__init__(self.path, cached_chunks=cached_chunks)
        try:
            if self.metadata.get("format_version") != "1":
                raise ValueError(
                    "Unsupported reference format; rebuild with 'fusion-function prepare-data'"
                )
            if self.metadata.get("species") != SPECIES:
                raise ValueError("Only human reference databases are supported")
            if self.metadata.get("assembly") != ASSEMBLY:
                raise ValueError("Only GRCh38 reference databases are supported")
            if release is not None and self.metadata.get("release") != str(release):
                raise ValueError(
                    f"Requested release {release}, but database has release {self.metadata.get('release')}"
                )
        except BaseException:
            self.close()
            raise


_local = threading.local()


def close_default_reference() -> None:
    """Release this thread's automatically cached connection and genome chunks."""
    reader = getattr(_local, "reader", None)
    if reader is not None:
        reader.close()
        del _local.reader
        del _local.signature


def get_reference(
    *,
    reference: ReferenceReader | None = None,
    database: str | Path | None = None,
    release: int | None = None,
) -> ReferenceReader:
    if reference is not None:
        if database is not None or release is not None:
            raise ValueError("Use reference, or database/release; these cannot be combined")
        return reference
    path = resolve_database(database, release=release)
    stat = path.stat()
    signature = (os.getpid(), str(path), stat.st_ino, stat.st_size, stat.st_mtime_ns)
    if getattr(_local, "signature", None) != signature:
        close_default_reference()
        _local.reader = ReferenceDatabase(path, release=release)
        _local.signature = signature
    elif release is not None and _local.reader.metadata.get("release") != str(release):
        raise ValueError("Requested release does not match the selected database")
    return _local.reader
