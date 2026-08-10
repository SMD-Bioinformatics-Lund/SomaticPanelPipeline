#!/usr/bin/env python3
"""Create a file manifest for a published Nextflow output directory."""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def is_excluded(relpath: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(relpath, pattern) for pattern in patterns)


def build_manifest(root: Path, case_id: str | None, excludes: list[str]) -> dict:
    root = root.resolve()
    if not root.is_dir():
        raise SystemExit(f"published directory does not exist: {root}")

    files = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relpath = path.relative_to(root).as_posix()
        if is_excluded(relpath, excludes):
            continue
        stat = path.stat()
        files.append(
            {
                "path": relpath,
                "size": stat.st_size,
                "sha256": sha256_file(path),
            }
        )

    return {
        "case_id": case_id,
        "root": str(root),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "file_count": len(files),
        "files": files,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-dir", required=True, help="Published output directory to scan.")
    parser.add_argument("--case-id", help="Case identifier from published_output_test_cases.json.")
    parser.add_argument("--out", required=True, help="Output manifest JSON.")
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="Relative-path glob to exclude. Repeatable.",
    )
    args = parser.parse_args()

    manifest = build_manifest(Path(args.published_dir), args.case_id, args.exclude)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
