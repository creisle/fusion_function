"""Command-line entry points."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from . import data


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "prepare-data":
        return data.main(arguments[1:])
    if arguments and arguments[0] == "export-reference":
        from . import prebuilt

        return prebuilt.main(arguments[1:])
    parser = argparse.ArgumentParser(prog="fusion-function")
    parser.add_argument(
        "command",
        choices=["prepare-data", "export-reference"],
        help="Install, build or export the human reference",
    )
    parser.parse_args(arguments)
    return 0
