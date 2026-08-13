#!/usr/bin/env python3
"""Compress stdin to a gzip file in a separate process."""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--compresslevel", type=int, choices=range(1, 10), default=1)
    args = parser.parse_args()
    with gzip.open(args.output, "wb", compresslevel=args.compresslevel) as output:
        shutil.copyfileobj(sys.stdin.buffer, output, length=1024 * 1024)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
