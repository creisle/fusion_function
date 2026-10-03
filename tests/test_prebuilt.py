"""Offline round trips and failure recovery for published reference artifacts."""

from __future__ import annotations

import copy
import gzip
import io
import json
import sqlite3
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pytest

from fusion_function import annotate_fusion_domains, data, prebuilt
from fusion_function.cli import main
from fusion_function.reference import ReferenceDatabase
from tests.preprocessing_fixture import fixture
from tests.test_prepare_options import entries_file
from tests import test_human_reference


class Response(io.BytesIO):
    def __init__(
        self, body: bytes, status: int = 200, headers: dict[str, str] | None = None
    ) -> None:
        super().__init__(body)
        self.status = status
        self.headers = headers or {}


@pytest.fixture
def exported(tmp_path: Path) -> tuple[Path, Path, dict[str, Any]]:
    source = tmp_path / "source.sqlite"
    db = fixture(source)
    db.execute("CREATE TABLE source_files (url TEXT PRIMARY KEY, sha256 TEXT, bytes INTEGER)")
    db.execute(
        "INSERT INTO source_files VALUES ('https://example.test/source.gz', 'source-hash', 123)"
    )
    data.preprocess_reference(db, entries_file(tmp_path))
    db.close()
    manifest = prebuilt.export_reference(source, tmp_path / "exports", zenodo_record=1234567)
    entry = json.loads(manifest.read_text())["references"][0]
    return source, manifest, entry


def mocked_download(
    monkeypatch: pytest.MonkeyPatch, manifest: Path, entry: dict[str, Any]
) -> list[str]:
    body = (manifest.parent / entry["filename"]).read_bytes()
    requests: list[str] = []

    def respond(request: urllib.request.Request) -> Response:
        requests.append(request.full_url)
        return Response(body)

    monkeypatch.setattr(prebuilt, "open_download", respond)
    return requests


def test_export_preserves_source_and_complete_annotation(
    exported: tuple[Path, Path, dict[str, Any]], tmp_path: Path
) -> None:
    source, manifest, entry = exported
    compact = tmp_path / "compact.sqlite"
    with gzip.open(manifest.parent / entry["filename"], "rb") as stream:
        compact.write_bytes(stream.read())
    assert prebuilt.data.sha256_file(compact) == entry["database_sha256"]
    with ReferenceDatabase(source) as original, ReferenceDatabase(compact) as copied:
        tables = {
            row[0] for row in copied.db.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }
        assert tables == set(prebuilt.RUNTIME_TABLES) | {"source_files"}
        assert (
            copied.db.execute("SELECT * FROM source_files").fetchall()
            == original.db.execute("SELECT * FROM source_files").fetchall()
        )
        assert original.db.execute("SELECT COUNT(*) FROM ensembl_transcript").fetchone()[0] == 4
        assert "reference_kind" not in original.metadata
        assert copied.metadata["reference_kind"] == "runtime"
        assert json.loads(copied.metadata["data_source_notices"]) == prebuilt.DATA_NOTICES
        for transcript_id in (
            "ENST00000000001",
            "ENST00000000002",
            "ENST00000000003",
            "ENST00000000004",
        ):
            assert original.get_transcript(transcript_id) == copied.get_transcript(transcript_id)
        arguments = dict(
            transcript1_id="ENST00000000001",
            transcript2_id="ENST00000000002",
            breakpoint1="1:111",
            breakpoint2="1:148",
            gene1_terminus="N",
            gene2_terminus="C",
        )
        assert annotate_fusion_domains(**arguments, reference=original) == annotate_fusion_domains(
            **arguments, reference=copied
        )
        assert original.sequence("pep", "ENSP00000000001") == copied.sequence(
            "pep", "ENSP00000000001"
        )


def test_human_model_checks_accept_runtime_export(
    exported: tuple[Path, Path, dict[str, Any]], tmp_path: Path
) -> None:
    """Functional built-reference checks must not require omitted source tables."""
    _, manifest, entry = exported
    compact = tmp_path / "compact.sqlite"
    with gzip.open(manifest.parent / entry["filename"], "rb") as stream:
        compact.write_bytes(stream.read())
    with sqlite3.connect(compact) as db:
        db.execute(
            "UPDATE ff_transcripts SET transcript_id='ENST00000305877' "
            "WHERE transcript_id='ENST00000000001'"
        )
    with ReferenceDatabase(compact) as reader:
        test_human_reference.test_real_transcript_exons_cds_and_peptide(
            reader, "BCR", "ENSG00000186716", "1", 1
        )
        test_human_reference.test_real_domains_have_interpro_metadata(
            reader, "BCR", "ENSG00000186716", "1", 1
        )
        with pytest.raises(pytest.skip.Exception, match="source coordinate tables"):
            test_human_reference.test_all_dna_records_match_grch38_core(reader)


def test_prepare_data_downloads_automatically_and_reuses_existing(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    _, manifest, entry = exported
    requests = mocked_download(monkeypatch, manifest, entry)
    monkeypatch.setattr(data, "open_url", lambda url: Response(manifest.read_bytes()))
    monkeypatch.setattr(
        data, "build", lambda args: pytest.fail("Source builder ran despite published reference")
    )
    cache = tmp_path / "cache"
    assert data.main(["--release", "116", "--cache-dir", str(cache)]) == 0
    output = cache / "homo_sapiens/GRCh38/release-116/ensembl.sqlite"
    with ReferenceDatabase(output) as reader:
        assert reader.metadata["reference_kind"] == "runtime"
    assert requests == [entry["url"]]
    assert data.main(["--release", "116", "--cache-dir", str(cache)]) == 0
    assert requests == [entry["url"]]
    assert not list(cache.rglob("*.part"))


@pytest.mark.parametrize("range_supported", [True, False])
def test_interrupted_download_resumes_or_restarts(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    range_supported: bool,
) -> None:
    _, manifest, entry = exported
    body = (manifest.parent / entry["filename"]).read_bytes()

    class Interrupted(Response):
        def read(self, size: int = -1) -> bytes:
            if self.tell():
                raise OSError("connection interrupted")
            return super().read(50)

    monkeypatch.setattr(prebuilt, "open_download", lambda request: Interrupted(body))
    cache = tmp_path / "cache"
    with pytest.raises(OSError, match="interrupted"):
        prebuilt.download_archive(entry, cache)
    partial = next(cache.rglob("*.part"))
    assert partial.read_bytes() == body[:50]

    def resumed(request: urllib.request.Request) -> Response:
        assert request.get_header("Range") == "bytes=50-"
        return (
            Response(body[50:], 206, {"Content-Range": f"bytes 50-{len(body) - 1}/{len(body)}"})
            if range_supported
            else Response(body)
        )

    monkeypatch.setattr(prebuilt, "open_download", resumed)
    assert prebuilt.download_archive(entry, cache).read_bytes() == body
    assert not partial.exists()


@pytest.mark.parametrize(
    "damage", ["archive_sha256", "database_sha256", "release", "database_bytes"]
)
def test_failed_install_preserves_previous_output(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    damage: str,
) -> None:
    _, manifest, original_entry = exported
    entry = copy.deepcopy(original_entry)
    if damage == "release":
        entry["release"] = 100  # Catalog claims a different release from the DB.
    elif damage.endswith("sha256"):
        entry[damage] = "0" * 64
    else:
        entry[damage] -= 1
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({"catalog_version": 1, "references": [entry]}))
    mocked_download(monkeypatch, manifest, original_entry)
    output = tmp_path / "existing.sqlite"
    output.write_bytes(b"previous reference")
    with pytest.raises(ValueError):
        prebuilt.install_reference(
            release=entry["release"],
            cache_dir=tmp_path / "cache",
            output=output,
            force=True,
            catalog=str(catalog),
        )
    assert output.read_bytes() == b"previous reference"
    assert not output.with_name(output.name + ".prebuilt.part").exists()


def test_latest_ensembl_is_not_replaced_with_older_prebuilt(
    exported: tuple[Path, Path, dict[str, Any]], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, manifest, entry = exported
    monkeypatch.setattr(data, "resolve_core", lambda base, release: (117, "", ""))
    monkeypatch.setattr(
        prebuilt, "open_download", lambda request: pytest.fail("Downloaded older release")
    )
    assert (
        prebuilt.install_reference(
            release=None, cache_dir=tmp_path / "cache", catalog=str(manifest)
        )
        is None
    )
    newer = copy.deepcopy(entry)
    newer["build_revision"] = 2
    assert prebuilt.select_reference([entry, newer], 116) == newer
    incompatible = copy.deepcopy(newer)
    incompatible["metadata"]["transcript_payload_codec"] = "unknown"
    assert prebuilt.select_reference([entry, incompatible], 116) == entry


def test_export_cli_and_unpublished_manifest(tmp_path: Path) -> None:
    source = tmp_path / "source.sqlite"
    db = fixture(source)
    data.preprocess_reference(db, entries_file(tmp_path))
    db.close()
    before = source.read_bytes()
    assert main(["export-reference", str(source), "--output", str(tmp_path / "exports")]) == 0
    manifest = next((tmp_path / "exports").glob("*.json"))
    assert prebuilt.read_catalog(str(manifest)) == []
    assert source.read_bytes() == before
    assert main(["export-reference", str(source), "--output", str(tmp_path / "exports")]) == 1


@pytest.mark.parametrize("release", [None, 100, 116])
def test_export_cli_selects_latest_or_pinned_local_release(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    release: int | None,
) -> None:
    """Omitting the path uses the configured cache and an optional local release."""
    source, _, _ = exported
    cache = tmp_path / "cache"
    monkeypatch.setenv("FUSION_FUNCTION_CACHEDIR", str(cache))
    monkeypatch.delenv("FUSION_FUNCTION_DB", raising=False)
    for number in (100, 116):
        target = cache / f"homo_sapiens/GRCh38/release-{number}/ensembl.sqlite"
        target.parent.mkdir(parents=True)
        target.write_bytes(source.read_bytes())
        with sqlite3.connect(target) as db:
            db.execute("UPDATE build_metadata SET value=? WHERE key='release'", (str(number),))
    output = tmp_path / "selected-export"
    arguments = ["export-reference", "--output", str(output), "--build-revision", "2"]
    if release is not None:
        arguments += ["--release", str(release)]
    assert main(arguments) == 0
    entry = json.loads(next(output.glob("*.json")).read_text())["references"][0]
    assert entry["release"] == (116 if release is None else release)
    assert entry["build_revision"] == 2


def test_export_cli_database_override_and_release_validation(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Honor the DB environment override; reject release mismatches and missing DBs."""
    source, _, _ = exported
    monkeypatch.setenv("FUSION_FUNCTION_DB", str(source))
    output = tmp_path / "override-export"
    assert main(["export-reference", "--output", str(output), "--release", "116"]) == 0
    assert main(["export-reference", "--output", str(output), "--release", "100"]) == 1
    assert "Requested release 100, but database has release 116" in caplog.text
    monkeypatch.setenv("FUSION_FUNCTION_DB", str(tmp_path / "missing.sqlite"))
    assert main(["export-reference", "--output", str(output)]) == 1
    assert "Reference database not found" in caplog.text
    # A positional path takes precedence over the environment setting.
    assert main(["export-reference", str(source), "--output", str(tmp_path / "explicit")]) == 0


def test_export_cli_rejects_invalid_release_before_selection(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["export-reference", "--output", str(tmp_path), "--release", "0"])
    assert exc.value.code == 2


def test_compact_reference_cannot_be_reprocessed(
    exported: tuple[Path, Path, dict[str, Any]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    _, manifest, entry = exported
    mocked_download(monkeypatch, manifest, entry)
    output = prebuilt.install_reference(
        release=116, cache_dir=tmp_path / "cache", catalog=str(manifest)
    )
    assert output is not None
    before = output.read_bytes()
    monkeypatch.setattr(
        data,
        "download",
        lambda url, directory: pytest.fail("Downloaded before validating compact DB"),
    )
    assert data.main(["--preprocess-only", str(output)]) == 1
    assert "--from-source" in caplog.text
    assert output.read_bytes() == before


def test_empty_catalog_and_explicit_source_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    catalog = tmp_path / "empty.json"
    catalog.write_text('{"catalog_version":1,"references":[]}')
    calls: list[bool] = []
    monkeypatch.setattr(data, "build", lambda args: calls.append(True) or tmp_path / "built.sqlite")
    assert (
        data.main(
            ["--reference-catalog", str(catalog), "--release", "116", "--cache-dir", str(tmp_path)]
        )
        == 0
    )
    monkeypatch.setattr(
        prebuilt,
        "read_catalog",
        lambda source=None: pytest.fail("Catalog read on explicit source build"),
    )
    assert data.main(["--from-source", "--cache-dir", str(tmp_path)]) == 0
    assert calls == [True, True]


def test_catalog_outage_uses_bundled_copy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def unavailable(url: str) -> Response:
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(data, "open_url", unavailable)
    assert prebuilt.read_catalog() == []
    with pytest.raises(urllib.error.URLError):
        prebuilt.read_catalog("https://example.test/catalog.json")


@pytest.mark.parametrize(
    "options",
    [
        ["--preprocess-only", "ref.sqlite", "--from-source"],
        ["--preprocess-only", "ref.sqlite", "--reference-catalog", "catalog.json"],
        ["--from-source", "--reference-catalog", "catalog.json"],
        ["--base-url", "https://mirror.test/pub/", "--reference-catalog", "catalog.json"],
    ],
)
def test_inapplicable_options_fail_before_network(options: list[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        data.main(options)
    assert exc.value.code == 2


@pytest.mark.parametrize(
    "field,value",
    [
        ("url", "http://example.test/db.gz"),
        ("url", "https://user:secret@example.test/db.gz"),
        ("archive_sha256", "invalid"),
        ("release", True),
        ("database_bytes", -1),
    ],
)
def test_invalid_catalog_entries_are_rejected(
    exported: tuple[Path, Path, dict[str, Any]], tmp_path: Path, field: str, value: object
) -> None:
    _, _, entry = exported
    entry[field] = value
    catalog = tmp_path / "bad.json"
    catalog.write_text(json.dumps({"catalog_version": 1, "references": [entry]}))
    with pytest.raises(ValueError):
        prebuilt.read_catalog(str(catalog))
