#!/usr/bin/env python3
"""Render baseline/candidate Nextflow command templates for parity cases."""

from __future__ import annotations

import argparse
import json
import os
import shlex
from pathlib import Path


def _quote(value: str) -> str:
    return shlex.quote(value)


def _case_command(case: dict, outdir_var: str, config: str, extra_args: list[str]) -> str:
    profiles = ",".join(case["profiles"])
    csv_value = "${" + case["csv_variable"] + "}"
    outdir_value = "${" + outdir_var + "}"
    parts = [
        "nextflow",
        "run",
        "main.nf",
        "-entry",
        case["entry"],
        "-profile",
        profiles,
        "-c",
        config,
        "--csv",
        csv_value,
        "--outdir",
        outdir_value,
        *extra_args,
    ]
    return " ".join(_quote(part) if not part.startswith("${") else part for part in parts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        default=".internal/published_output_test_cases.json",
        help="JSON test-case manifest.",
    )
    parser.add_argument(
        "--config",
        default="configs/nextflow.hopper.config",
        help="Nextflow config path to render into commands.",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Extra argument to append to every rendered Nextflow command. Repeatable.",
    )
    args = parser.parse_args()

    with open(args.cases) as fh:
        manifest = json.load(fh)

    print("# Published Output Parity Commands")
    print()
    print("Set the CSV and output variables described in the manifest before running these commands.")
    print("Run baseline commands on main/master and candidate commands on the containers branch.")
    print()
    for case in manifest["cases"]:
        print(f"## {case['id']}")
        print(f"# {case['reason']}")
        print("# baseline")
        print(_case_command(case, "BASELINE_OUTDIR", args.config, args.extra_arg))
        print("# candidate")
        print(_case_command(case, "CANDIDATE_OUTDIR", args.config, args.extra_arg))
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
