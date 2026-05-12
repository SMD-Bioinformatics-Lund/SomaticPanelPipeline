#!/usr/bin/env python3
"""Filter paired TNscope VCF records for somatic variant review."""

from argparse import ArgumentParser

import vcf2 as vcf2
from vcf_pipeline_utils import (  # noqa: E402
    append_filter_status,
    leading_float,
    open_output,
    write_filter_header,
    write_variant,
)

MIN_VAF_RATIO = 3
MIN_VAF_HOMOPOLYMER_RATIO = 5
FILTER_HEADERS = [
    '##FILTER=<ID=FAIL_NVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_LOW_TCOV,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_VERYLOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_LOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_INDEL,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_SNV,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_HOMOPOLYMER_INDEL,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_HOMOPOLYMER_SNV,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_STRANDBIAS,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_STRANDBIAS,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_PVALUE,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_NO_TVAR,Description="Record fails the filters">',
    '##FILTER=<ID=.,Description="Record fails the filters">',
]


def first_two_counts(value):
    """Return the first two comma-delimited count values.

    Args:
        value: Allelic depth string.

    Returns:
        Tuple of reference and alternate depths.
    """
    parts = str(value or "").split(",")
    ref_count = leading_float(parts[0]) if len(parts) > 0 else 0.0
    alt_count = leading_float(parts[1]) if len(parts) > 1 else 0.0
    return ref_count, alt_count


def filter_vcf(
    vcf_file,
    tumor_sample,
    normal_sample,
    out_file,
    min_vaf_ratio=MIN_VAF_RATIO,
    min_vaf_homopolymer_ratio=MIN_VAF_HOMOPOLYMER_RATIO,
    min_tumor_depth=100,
    very_low_tumor_vaf=0.02,
    low_tumor_vaf=0.05,
    max_homopolymer_unit_length=2,
    min_homopolymer_repeat_count=10,
    indel_warn_sor=4,
    indel_fail_sor=10,
    snv_warn_sor=2.5,
    snv_fail_sor=4,
    max_pvalue=0.05,
):
    """Filter paired TNscope VCF records.

    Args:
        vcf_file: Input VCF path.
        tumor_sample: Tumor sample ID.
        normal_sample: Normal sample ID.
        out_file: Output VCF path.
    """
    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        write_filter_header(reader, out_fh, FILTER_HEADERS)
        for var in reader:
            statuses = []
            vaf = {}
            depth = {}
            calculated_vaf = {}
            is_indel = len(var["REF"]) != len(var["ALT"])

            repeat_units_num = 0.0
            repeat_seq = ""
            if var["INFO"].get("RPA") and var["INFO"].get("RU"):
                repeat_units_num = min(leading_float(item) for item in var["INFO"]["RPA"].split(","))
                repeat_seq = var["INFO"].get("RU") or ""

            for gt in var["GT"]:
                sample_type = "N" if gt["_sample_id"] == normal_sample else "T"
                vaf[sample_type] = leading_float(gt.get("AF"))
                ref_count, alt_count = first_two_counts(gt.get("AD"))
                sample_depth = ref_count + alt_count
                if sample_depth > 0:
                    depth[sample_type] = sample_depth
                    calculated_vaf[sample_type] = float(f"{alt_count / sample_depth:.3f}")
                else:
                    depth[sample_type] = 0.0
                    calculated_vaf[sample_type] = 0.0

            if vaf.get("T", 0) > 0:
                if not calculated_vaf.get("N"):
                    calculated_vaf["N"] = 1e-10
                if calculated_vaf["N"] > 0 and calculated_vaf.get("T", 0) / calculated_vaf["N"] < min_vaf_ratio:
                    statuses.append("FAIL_NVAF")
                if depth.get("T", 0) < min_tumor_depth:
                    statuses.append("WARN_LOW_TCOV")
                if vaf["T"] < very_low_tumor_vaf:
                    statuses.append("WARN_VERYLOW_TVAF")
                elif vaf["T"] < low_tumor_vaf:
                    statuses.append("WARN_LOW_TVAF")
                if len(repeat_seq) <= max_homopolymer_unit_length and repeat_units_num >= min_homopolymer_repeat_count:
                    if is_indel or not vaf.get("N", 0) or vaf["T"] / vaf["N"] >= min_vaf_homopolymer_ratio:
                        statuses.append("WARN_HOMOPOLYMER_" + ("INDEL" if is_indel else "SNV"))
                    else:
                        statuses.append("FAIL_HOMOPOLYMER_" + ("INDEL" if is_indel else "SNV"))

                sor = leading_float(var["INFO"].get("SOR"))
                if is_indel and sor > indel_warn_sor:
                    statuses.append("FAIL_STRANDBIAS" if sor > indel_fail_sor else "WARN_STRANDBIAS")
                if not is_indel and sor > snv_warn_sor:
                    statuses.append("FAIL_STRANDBIAS" if sor > snv_fail_sor else "WARN_STRANDBIAS")
                if leading_float(var["INFO"].get("PV2")) > max_pvalue:
                    statuses.append("FAIL_PVALUE")
            else:
                statuses.append("FAIL_NO_TVAR")

            append_filter_status(var, statuses)
            write_variant(var, out_fh)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Filter paired TNscope VCF records.")
    parser.add_argument("--vcf", required=True, help="Input TNscope VCF file.")
    parser.add_argument("--tumor", required=True, help="Tumor sample name.")
    parser.add_argument("--normal", required=True, help="Normal sample name.")
    parser.add_argument("--out", required=True, help="Output filtered VCF file.")
    add_filter_threshold_args(parser)
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    filter_vcf(
        args.vcf,
        args.tumor,
        args.normal,
        args.out,
        args.min_vaf_ratio,
        args.min_vaf_homopolymer_ratio,
        args.min_tumor_depth,
        args.very_low_tumor_vaf,
        args.low_tumor_vaf,
        args.max_homopolymer_unit_length,
        args.min_homopolymer_repeat_count,
        args.indel_warn_sor,
        args.indel_fail_sor,
        args.snv_warn_sor,
        args.snv_fail_sor,
        args.max_pvalue,
    )
    return 0


def add_filter_threshold_args(parser):
    """Add configurable TNscope paired filtering thresholds."""
    parser.add_argument("--min-vaf-ratio", type=float, default=MIN_VAF_RATIO)
    parser.add_argument("--min-vaf-homopolymer-ratio", type=float, default=MIN_VAF_HOMOPOLYMER_RATIO)
    parser.add_argument("--min-tumor-depth", type=float, default=100)
    parser.add_argument("--very-low-tumor-vaf", type=float, default=0.02)
    parser.add_argument("--low-tumor-vaf", type=float, default=0.05)
    parser.add_argument("--max-homopolymer-unit-length", type=float, default=2)
    parser.add_argument("--min-homopolymer-repeat-count", type=float, default=10)
    parser.add_argument("--indel-warn-sor", type=float, default=4)
    parser.add_argument("--indel-fail-sor", type=float, default=10)
    parser.add_argument("--snv-warn-sor", type=float, default=2.5)
    parser.add_argument("--snv-fail-sor", type=float, default=4)
    parser.add_argument("--max-pvalue", type=float, default=0.05)


if __name__ == "__main__":
    raise SystemExit(main())
