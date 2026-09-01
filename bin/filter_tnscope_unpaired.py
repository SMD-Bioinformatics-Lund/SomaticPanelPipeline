#!/usr/bin/env python3
"""Filter unpaired TNscope VCF records for variant review."""

from argparse import ArgumentParser

import vcf2 as vcf2
from vcf_pipeline_utils import (  # noqa: E402
    append_filter_status,
    leading_float,
    open_output,
    write_filter_header,
    write_variant,
)

FILTER_HEADERS = [
    '##FILTER=<ID=WARN_LOW_TCOV,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_VERYLOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_LOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_INDEL,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_SNV,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_STRANDBIAS,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_STRANDBIAS,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_PVALUE,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_NO_TVAR,Description="Record fails the filters">',
    '##FILTER=<ID=.,Description="Record fails the filters">',
]


def sum_counts(value):
    """Sum comma-delimited numeric count values.

    Args:
        value: Allelic depth string.

    Returns:
        Total count.
    """
    return sum(leading_float(item) for item in str(value or "").split(","))


def filter_vcf(
    vcf_file,
    out_file,
    min_tumor_depth=100,
    very_low_tumor_vaf=0.02,
    low_tumor_vaf=0.05,
    max_homopolymer_unit_length=2,
    min_homopolymer_repeat_count=10,
    indel_warn_sor=4,
    indel_fail_sor=10,
    snv_warn_sor=2.5,
    snv_fail_sor=4,
):
    """Filter unpaired TNscope VCF records.

    Args:
        vcf_file: Input VCF path.
        out_file: Output VCF path.
    """
    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        write_filter_header(reader, out_fh, FILTER_HEADERS)
        for var in reader:
            statuses = []
            vaf = {}
            depth = {}
            is_indel = len(var["REF"]) != len(var["ALT"])

            repeat_units_num = 0.0
            repeat_seq = ""
            if var["INFO"].get("RPA") and var["INFO"].get("RU"):
                repeat_units_num = min(leading_float(item) for item in var["INFO"]["RPA"].split(","))
                repeat_seq = var["INFO"].get("RU") or ""

            for gt in var["GT"]:
                vaf["T"] = leading_float(gt.get("AF"))
                depth["T"] = sum_counts(gt.get("AD"))

            if vaf.get("T", 0) > 0:
                if depth.get("T", 0) < min_tumor_depth:
                    statuses.append("WARN_LOW_TCOV")
                if vaf["T"] < very_low_tumor_vaf:
                    statuses.append("WARN_VERYLOW_TVAF")
                elif vaf["T"] < low_tumor_vaf:
                    statuses.append("WARN_LOW_TVAF")
                if len(repeat_seq) <= max_homopolymer_unit_length and repeat_units_num >= min_homopolymer_repeat_count:
                    statuses.append("WARN_HOMOPOLYMER_" + ("INDEL" if is_indel else "SNV"))

                sor = leading_float(var["INFO"].get("SOR"))
                if is_indel and sor > indel_warn_sor:
                    statuses.append("FAIL_STRANDBIAS" if sor > indel_fail_sor else "WARN_STRANDBIAS")
                if not is_indel and sor > snv_warn_sor:
                    statuses.append("FAIL_STRANDBIAS" if sor > snv_fail_sor else "WARN_STRANDBIAS")
            else:
                statuses.append("FAIL_NO_TVAR")

            append_filter_status(var, statuses)
            write_variant(var, out_fh)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Filter unpaired TNscope VCF records.")
    parser.add_argument("--vcf", required=True, help="Input TNscope VCF file.")
    parser.add_argument("--out", required=True, help="Output filtered VCF file.")
    parser.add_argument("--min-tumor-depth", type=float, default=100)
    parser.add_argument("--very-low-tumor-vaf", type=float, default=0.02)
    parser.add_argument("--low-tumor-vaf", type=float, default=0.05)
    parser.add_argument("--max-homopolymer-unit-length", type=float, default=2)
    parser.add_argument("--min-homopolymer-repeat-count", type=float, default=10)
    parser.add_argument("--indel-warn-sor", type=float, default=4)
    parser.add_argument("--indel-fail-sor", type=float, default=10)
    parser.add_argument("--snv-warn-sor", type=float, default=2.5)
    parser.add_argument("--snv-fail-sor", type=float, default=4)
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    filter_vcf(
        args.vcf,
        args.out,
        args.min_tumor_depth,
        args.very_low_tumor_vaf,
        args.low_tumor_vaf,
        args.max_homopolymer_unit_length,
        args.min_homopolymer_repeat_count,
        args.indel_warn_sor,
        args.indel_fail_sor,
        args.snv_warn_sor,
        args.snv_fail_sor,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
