"""Create SNV panel-of-normals summaries from caller VCF files."""

from __future__ import annotations

from argparse import ArgumentParser
from glob import glob
from math import sqrt

import vcf2 as vcf
from vcf_pipeline_utils import leading_float, open_output
from aggregate_vcf import which_variantcaller


def average(values):
    """Return the arithmetic mean of values."""
    return sum(values) / len(values) if values else 0


def std_dev(values):
    """Return the population standard deviation using the legacy argument order."""
    if not values:
        return 0
    avg = values[0]
    return sqrt(sum((value - avg) ** 2 for value in values[1:]) / len(values))


def count_germline(values):
    """Count values that look germline in the PON input set."""
    return sum(1 for value in values if (0.35 < value < 0.65) or value > 0.85)


def mean_nongermline(values):
    """Summarize non-germline VAF values."""
    kept = [value for value in values if value < 0.35]
    return average(kept), std_dev(kept), ",".join(str(value) for value in kept)


def first_two_ad(value):
    """Return ref and alt depths from an AD field."""
    fields = str(value or "").split(",")
    ref = leading_float(fields[0]) if fields else 0
    alt = leading_float(fields[1]) if len(fields) > 1 else 0
    return ref, alt


def get_stats(var, caller):
    """Return VAF, variant depth, and total depth for the first sample."""
    if caller in {"mutect2", "tnscope", "vardict", "pindel"}:
        for gt in var["GT"]:
            ref_depth, alt_depth = first_two_ad(gt.get("AD"))
            depth = ref_depth + alt_depth
            af = alt_depth / depth if depth > 0 else 0
            return leading_float(gt.get("AF")) or af, alt_depth, depth
    if caller == "freebayes":
        for gt in var["GT"]:
            vaf = 0
            vd = 0
            if gt.get("AO") and gt.get("AO") != ".":
                vaf = float(f"{leading_float(gt.get('AO')) / leading_float(gt.get('DP')):.4f}")
                vd = leading_float(gt.get("AO"))
            return vaf, vd, leading_float(gt.get("DP"))
    return 0, 0, 0


def create_pon(filemask, out_file):
    """Create a PON summary table.

    Args:
        filemask: Glob pattern for input VCF files.
        out_file: Output table path, or ``-`` for stdout.
    """
    files = glob(filemask)
    all_values = {}
    for path in files:
        with vcf.VCFReader(path) as reader:
            caller = which_variantcaller(reader.meta)
            for var in reader:
                vaf_value, _vd, depth = get_stats(var, caller)
                if depth == 0:
                    continue
                varid = f"{var['CHROM']}_{var['POS']}_{var['REF']}_{var['ALT']}"
                all_values.setdefault(varid, []).append(vaf_value)

    with open_output(out_file) as out_fh:
        for varid in sorted(all_values):
            values = all_values[varid]
            num_vars = len(values)
            num_germline = count_germline(values)
            mean, sd, kept_values = mean_nongermline(values)
            if num_germline == num_vars:
                continue
            out_fh.write(f"{varid}\t{num_vars}\t{num_germline}\t{len(files)}\t{mean}\t{sd}\t{kept_values}\n")


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Create an SNV panel-of-normals table.")
    parser.add_argument("--vcf-mask", required=True, help="Input VCF glob pattern.")
    parser.add_argument("--out", required=True, help="Output PON table.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    create_pon(args.vcf_mask, args.out)
    return 0
