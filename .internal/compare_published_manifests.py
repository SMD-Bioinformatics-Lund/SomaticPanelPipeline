#!/usr/bin/env python3
"""Compare two published-output manifests."""

from __future__ import annotations

import argparse
import fnmatch
import json
from pathlib import Path


def load_manifest(path: Path) -> dict:
    with path.open() as fh:
        manifest = json.load(fh)
    return {item["path"]: item for item in manifest["files"]}


def ignored(path: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(path, pattern) for pattern in patterns)


def compare(baseline: dict, candidate: dict, ignore: list[str]) -> dict:
    base_paths = {path for path in baseline if not ignored(path, ignore)}
    cand_paths = {path for path in candidate if not ignored(path, ignore)}
    common = sorted(base_paths & cand_paths)

    changed = []
    unchanged = []
    for path in common:
        base = baseline[path]
        cand = candidate[path]
        if base["sha256"] == cand["sha256"] and base["size"] == cand["size"]:
            unchanged.append(path)
        else:
            changed.append(
                {
                    "path": path,
                    "baseline_size": base["size"],
                    "candidate_size": cand["size"],
                    "baseline_sha256": base["sha256"],
                    "candidate_sha256": cand["sha256"],
                }
            )

    return {
        "missing_in_candidate": sorted(base_paths - cand_paths),
        "extra_in_candidate": sorted(cand_paths - base_paths),
        "changed": changed,
        "unchanged": unchanged,
    }


def write_markdown(report: dict, out_path: Path) -> None:
    lines = [
        "# Published Output Comparison",
        "",
        f"- Missing in candidate: {len(report['missing_in_candidate'])}",
        f"- Extra in candidate: {len(report['extra_in_candidate'])}",
        f"- Changed content: {len(report['changed'])}",
        f"- Unchanged: {len(report['unchanged'])}",
        "",
    ]

    sections = [
        ("Missing In Candidate", report["missing_in_candidate"]),
        ("Extra In Candidate", report["extra_in_candidate"]),
        ("Changed Content", report["changed"]),
    ]
    for title, rows in sections:
        lines.extend([f"## {title}", ""])
        if not rows:
            lines.extend(["None.", ""])
            continue
        if title == "Changed Content":
            lines.extend(["| Path | Baseline Size | Candidate Size |", "|---|---:|---:|"])
            for row in rows:
                lines.append(f"| `{row['path']}` | {row['baseline_size']} | {row['candidate_size']} |")
        else:
            for path in rows:
                lines.append(f"- `{path}`")
        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, help="Baseline manifest JSON.")
    parser.add_argument("--candidate", required=True, help="Candidate manifest JSON.")
    parser.add_argument("--json-out", required=True, help="Output comparison JSON.")
    parser.add_argument("--md-out", required=True, help="Output comparison Markdown.")
    parser.add_argument(
        "--ignore",
        action="append",
        default=[],
        help="Relative-path glob to ignore during comparison. Repeatable.",
    )
    args = parser.parse_args()

    report = compare(load_manifest(Path(args.baseline)), load_manifest(Path(args.candidate)), args.ignore)

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    with json_out.open("w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
        fh.write("\n")
    write_markdown(report, Path(args.md_out))

    return 1 if report["missing_in_candidate"] or report["extra_in_candidate"] or report["changed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
