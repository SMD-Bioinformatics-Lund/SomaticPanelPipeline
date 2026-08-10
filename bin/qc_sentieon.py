#!/usr/bin/env python3
"""Parse Sentieon QC metric files in the current directory."""

import json
from argparse import ArgumentParser

from vcf_pipeline_utils import leading_float, open_output


def qc_sentieon(sample_id, assay_type, out_file):
    """Parse Sentieon QC metric files in the current directory."""
    results = {}
    pct_above = {}
    if assay_type == "panel":
        pct_above, median = panel_coverage_calc()
        results["median_cov"] = median
        for row in metric_value_rows("hs_metrics.txt"):
            results["pct_on_target"] = row[18]
            results["fold_enrichment"] = row[25]
            results["mean_coverage"] = row[22]
            results["fold_80"] = row[32]
            results["at_drop"] = row[48]
            results["gc_drop"] = row[49]
    elif assay_type == "wgs":
        pct_above = parse_wgs_metrics(results)
        parse_gc_summary(results)
    else:
        raise ValueError("assay type must be panel or wgs")

    for row in metric_value_rows("is_metrics.txt"):
        results["ins_size"] = f"{leading_float(row[4]):.0f}"
        results["ins_size_dev"] = f"{leading_float(row[5]):.0f}"
    for row in metric_value_rows("dedup_metrics.txt"):
        results["dup_reads"] = row[6]
        results["num_reads"] = row[2]
        results["dup_pct"] = row[8]
        results["mapped_reads"] = leading_float(row[2]) - leading_float(row[4])
    for row in metric_value_rows("aln_metrics.txt"):
        results["pf_mismatch_rate"] = row[12]
        results["pf_error_rate"] = row[13]

    results["pct_above_x"] = pct_above
    results["sample_id"] = sample_id
    with open_output(out_file) as out_fh:
        json.dump(results, out_fh, indent=3)
        out_fh.write("\n")


def metric_value_rows(path):
    """Yield Sentieon metric value rows following command/header lines."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#SentieonCommandLine"):
                next(fh, None)
                yield next(fh).rstrip("\n").split("\t")


def panel_coverage_calc():
    """Calculate panel coverage summary values."""
    cov = []
    with open("cov_metrics.txt") as fh:
        for line in fh:
            if not line.startswith("Locus"):
                vals = next(fh, "").split("\t")
                if len(vals) > 1:
                    cov.append(leading_float(vals[1]))

    sorted_cov = sorted(cov)
    mid = len(sorted_cov) // 2
    median = sorted_cov[mid] if len(sorted_cov) % 2 else sorted_cov[mid - 1]

    pct = {}
    with open("cov_metrics.txt.sample_summary") as fh:
        for line in fh:
            if line.startswith("sample_id"):
                vals = next(fh).rstrip("\n").split("\t")
                for key, idx in (("1", 6), ("10", 7), ("30", 8), ("100", 9), ("250", 10), ("500", 11)):
                    pct[key] = vals[idx]
    return pct, median


def parse_wgs_metrics(results):
    """Parse WGS metrics and return pct-above thresholds."""
    pct = {}
    quartiles = {}
    total = 0
    pct_obs = {}
    with open("wgs_metrics.txt") as fh:
        for line in fh:
            if line.startswith("#SentieonCommandLine"):
                next(fh, None)
                row = next(fh).split("\t")
                results["median_cov"] = row[3]
                results["sd_coverage"] = row[2]
                results["mean_coverage"] = row[1]
                pct_obs = {
                    "R_25": (leading_float(row[0]) + 1) / 4,
                    "R_50": (leading_float(row[0]) + 1) / 2,
                    "R_75": ((leading_float(row[0]) + 1) * 3) / 4,
                }
                for key, idx in zip(
                    ("1", "5", "10", "15", "20", "25", "30", "40", "50", "60", "70", "80", "90", "100"),
                    range(12, 26),
                ):
                    pct[key] = row[idx]
            elif line[:1].isdigit():
                depth, count = line.split()[:2]
                total += leading_float(count)
                for key, obs in pct_obs.items():
                    if total >= obs and key not in quartiles:
                        quartiles[key] = leading_float(depth)
    results["iqr"] = quartiles.get("R_75", 0) - quartiles.get("R_25", 0)
    return pct


def parse_gc_summary(results):
    """Parse GC summary file into results."""
    for row in metric_value_rows("gc_summary.txt"):
        results["at_drop"] = row[5]
        results["gc_drop"] = row[6]


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--type", required=True, choices=["panel", "wgs"], help="Assay type.")
    parser.add_argument("--out", required=True, help="Output JSON.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    qc_sentieon(args.sample_id, args.type, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
