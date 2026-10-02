import os
from pathlib import Path

import pytest

from fusion_function.reference import ReferenceDatabase, close_default_reference, resolve_database


def pytest_addoption(parser):
    parser.addoption(
        "--human-reference-db",
        metavar="SQLITE",
        default=None,
        help="Override the automatically discovered human GRCh38 reference database",
    )


def pytest_generate_tests(metafunc):
    """Run original integration cases against both cached responses and a build."""
    if (
        "reference_db" in metafunc.fixturenames
        and metafunc.definition.get_closest_marker("integration") is not None
    ):
        metafunc.parametrize(
            "reference_db",
            [
                pytest.param("bundled", id="bundled"),
                pytest.param("built", id="built", marks=pytest.mark.human_reference),
            ],
            indirect=True,
            scope="session",
        )


@pytest.fixture(scope="session")
def human_reference_db(request):
    """Use an existing local build; never download or build a reference for tests."""
    supplied = request.config.getoption("--human-reference-db") or os.environ.get(
        "FUSION_FUNCTION_DB"
    )
    if supplied:
        path = Path(supplied).expanduser().resolve()
    else:
        try:
            path = resolve_database()
        except FileNotFoundError:
            pytest.skip(
                "No built human GRCh38 reference in the default cache; "
                "run prepare-data or supply --human-reference-db PATH"
            )
    if not path.is_file():
        pytest.fail(f"Human reference database does not exist: {path}")
    with ReferenceDatabase(path) as reference:
        # A sparse regression fixture must never pass as a full human reference.
        ready = reference.db.execute(
            "SELECT COUNT(*) FROM ff_transcripts WHERE status='ready'"
        ).fetchone()[0]
        if ready < 10000:
            pytest.fail(
                f"Expected a fully built human reference, found only {ready} ready transcripts in {path}"
            )
        yield reference


@pytest.fixture(scope="session")
def bundled_reference_db():
    """Keep the original response-derived reference available without a build."""
    with ReferenceDatabase(Path(__file__).parent / "data/reference.sqlite") as reference:
        yield reference


@pytest.fixture(scope="session")
def reference_db(request):
    """Select the parameterized reference; lightweight tests use the bundle."""
    fixture = (
        "human_reference_db"
        if getattr(request, "param", "bundled") == "built"
        else "bundled_reference_db"
    )
    return request.getfixturevalue(fixture)


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    import socket
    import urllib.request

    def blocked_connection(
        address, timeout=socket._GLOBAL_DEFAULT_TIMEOUT, source_address=None, *, all_errors=False
    ):
        raise AssertionError("Runtime/test network access is forbidden")

    def blocked_urlopen(
        url,
        data=None,
        timeout=socket._GLOBAL_DEFAULT_TIMEOUT,
        *,
        cafile=None,
        capath=None,
        cadefault=False,
        context=None,
    ):
        raise AssertionError("Runtime/test network access is forbidden")

    monkeypatch.setattr(socket, "create_connection", blocked_connection)
    monkeypatch.setattr(urllib.request, "urlopen", blocked_urlopen)
    yield
    close_default_reference()
