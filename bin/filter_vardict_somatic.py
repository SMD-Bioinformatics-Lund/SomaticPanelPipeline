#!/usr/bin/env python3
"""Filter paired VarDict VCF records for somatic variant review."""

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
    '##FILTER=<ID=FAIL_LONGDEL,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_LONGINS,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_LOW_TCOV,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_VERYLOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_LOW_TVAF,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_STRANDBIAS,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_MQ,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_INDEL,Description="Record fails the filters">',
    '##FILTER=<ID=WARN_HOMOPOLYMER_SNV,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_HOMOPOLYMER_INDEL,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_HOMOPOLYMER_SNV,Description="Record fails the filters">',
    '##FILTER=<ID=FAIL_NO_TVAR,Description="Record fails the filters">',
    '##FILTER=<ID=.,Description="Record fails the filters">',
]


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
    min_strand_bias=0.05,
    min_mapping_quality=10,
    max_homopolymer_unit_length=2,
    min_homopolymer_repeat_count=10,
    max_indel_length=250,
):
    """Filter paired VarDict VCF records.

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
            strand_bias = {}
            mapping_qual = {}

            msi = leading_float(var["INFO"].get("MSI"))
            msilen = leading_float(var["INFO"].get("MSILEN"))
            is_indel = len(var["REF"]) != len(var["ALT"])

            for gt in var["GT"]:
                sample_type = "N" if gt["_sample_id"] == normal_sample else "T"
                vaf[sample_type] = leading_float(gt.get("AF"))
                depth[sample_type] = leading_float(gt.get("DP"))
                strand_bias[sample_type] = leading_float(gt.get("SBF"))
                mapping_qual[sample_type] = leading_float(gt.get("MQ"))

            if len(var["REF"]) - len(var["ALT"]) > max_indel_length or var["ALT"] == "<DEL>":
                statuses.append("FAIL_LONGDEL")
            elif len(var["ALT"]) - len(var["REF"]) > max_indel_length or var["ALT"] in ("<DUP>", "<INS>"):
                statuses.append("FAIL_LONGINS")

            if vaf.get("T", 0) > 0:
                if vaf.get("N", 0) > 0 and vaf["T"] / vaf["N"] < min_vaf_ratio:
                    statuses.append("FAIL_NVAF")
                if depth.get("T", 0) < min_tumor_depth:
                    statuses.append("WARN_LOW_TCOV")
                if vaf["T"] < very_low_tumor_vaf:
                    statuses.append("WARN_VERYLOW_TVAF")
                elif vaf["T"] < low_tumor_vaf:
                    statuses.append("WARN_LOW_TVAF")
                if strand_bias.get("T", 0) < min_strand_bias:
                    statuses.append("WARN_STRANDBIAS")
                if mapping_qual.get("T", 0) <= min_mapping_quality:
                    statuses.append("WARN_MQ")
                if msilen <= max_homopolymer_unit_length and msi > min_homopolymer_repeat_count:
                    if not is_indel or not vaf.get("N", 0) or vaf["T"] / vaf["N"] >= min_vaf_homopolymer_ratio:
                        statuses.append("WARN_HOMOPOLYMER_" + ("INDEL" if is_indel else "SNV"))
                    else:
                        statuses.append("FAIL_HOMOPOLYMER_" + ("INDEL" if is_indel else "SNV"))
            else:
                statuses.append("FAIL_NO_TVAR")

            append_filter_status(var, statuses)
            write_variant(var, out_fh)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Filter paired VarDict VCF records.")
    parser.add_argument("--vcf", required=True, help="Input VarDict VCF file.")
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
        args.min_strand_bias,
        args.min_mapping_quality,
        args.max_homopolymer_unit_length,
        args.min_homopolymer_repeat_count,
        args.max_indel_length,
    )
    return 0


def add_filter_threshold_args(parser):
    """Add configurable VarDict paired filtering thresholds."""
    parser.add_argument("--min-vaf-ratio", type=float, default=MIN_VAF_RATIO)
    parser.add_argument("--min-vaf-homopolymer-ratio", type=float, default=MIN_VAF_HOMOPOLYMER_RATIO)
    parser.add_argument("--min-tumor-depth", type=float, default=100)
    parser.add_argument("--very-low-tumor-vaf", type=float, default=0.02)
    parser.add_argument("--low-tumor-vaf", type=float, default=0.05)
    parser.add_argument("--min-strand-bias", type=float, default=0.05)
    parser.add_argument("--min-mapping-quality", type=float, default=10)
    parser.add_argument("--max-homopolymer-unit-length", type=float, default=2)
    parser.add_argument("--min-homopolymer-repeat-count", type=float, default=10)
    parser.add_argument("--max-indel-length", type=float, default=250)


if __name__ == "__main__":
    raise SystemExit(main())
