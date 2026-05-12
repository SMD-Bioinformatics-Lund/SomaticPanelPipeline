#!/usr/bin/env python3
"""Filter paired FreeBayes VCF records for somatic variant review.

The script reads a VCF with tumor and normal samples, appends FreeBayes-specific
filter tags, adds an SSC INFO value, and writes the filtered VCF to a file.
"""

import re
import sys
from argparse import ArgumentParser
from contextlib import nullcontext

import vcf2 as vcf2

SSC_THRES = 30
LOD_THRES = 10
MIN_VAF_RATIO = 3
LEADING_FLOAT_RE = re.compile(r"^[\s]*([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)")


def leading_float(value):
    """Convert a scalar-like value to a number using leading numeric text.

    Args:
        value: Input value from a VCF field.

    Returns:
        Floating-point numeric value. Non-numeric or empty values become 0.0.
    """
    if value is None or value == "":
        return 0.0
    match = LEADING_FLOAT_RE.match(str(value))
    if match:
        return float(match.group(1))
    return 0.0


def format_compact_number(value):
    """Format a numeric value for VCF INFO output.

    Args:
        value: Numeric value to format.

    Returns:
        Compact decimal string without unnecessary trailing zeros.
    """
    return f"{value:.15g}"


def add_info(var, key, val):
    """Append an INFO field to a parsed VCF record.

    Args:
        var: Parsed variant dictionary.
        key: INFO key to append.
        val: INFO value to store.
    """
    var["INFO_order"].append(key)
    var["INFO"][key] = val


def open_output(path):
    """Open the output VCF destination.

    Args:
        path: Output file path, or ``-`` to write to stdout.

    Returns:
        Context manager yielding a writable text file handle.
    """
    if path == "-":
        return nullcontext(sys.stdout)
    return open(path, "w")


def vcfstr(var, out_fh):
    """Write a parsed VCF record.

    Args:
        var: Parsed variant dictionary to serialize.
        out_fh: Writable output file handle.
    """
    fields = [
        var["CHROM"],
        var["POS"],
        var["ID"],
        var["REF"],
        var["ALT"],
        var["QUAL"],
        var["FILTER"],
        ";".join(f"{key}={var['INFO'][key]}" for key in var["INFO_order"]),
        ":".join(var["FORMAT"]),
    ]
    fields.extend(":".join(gt.get(key, "") for key in var["FORMAT"]) for gt in var["GT"])
    out_fh.write("\t".join(fields) + "\t\n")


def print_header(reader, out_fh):
    """Write the original VCF header plus somatic filter metadata.

    Args:
        reader: Open VCF reader with parsed header lines.
        out_fh: Writable output file handle.
    """
    for line in reader.meta_header_lines:
        print(line, file=out_fh)

    print('##INFO=<ID=SSC,Number=1,Type=Float,Description="Somatic score">', file=out_fh)
    print('##FILTER=<ID=FAIL_GT,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=LOH,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=FAIL_QUAL,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=FAIL_LOD,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=WARN_NOVAR,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=FAIL_NVAF,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=WARN_MQ,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=.,Description="Record fails the filters">', file=out_fh)
    if reader.column_header_line:
        print(reader.column_header_line, file=out_fh)


def filter_vcf(
    vcf_file,
    tumor_sample,
    normal_sample,
    out_file,
    ssc_threshold=SSC_THRES,
    lod_threshold=LOD_THRES,
    min_vaf_ratio=MIN_VAF_RATIO,
    min_mapping_quality=10,
):
    """Filter a paired tumor/normal FreeBayes VCF.

    Args:
        vcf_file: Path to the input VCF file.
        tumor_sample: Tumor sample name from the VCF header.
        normal_sample: Normal sample name from the VCF header.
        out_file: Path to the output VCF file.
    """
    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        print_header(reader, out_fh)
        for var in reader:
            likelihood = {}
            gl_idx = {}
            genotype = {}
            altobs = {}
            depth = {}

            status = "PASS"
            for gt in var["GT"]:
                sample_type = "T"
                if gt["_sample_id"] == normal_sample:
                    sample_type = "N"

                if sample_type == "T" and gt.get("GT") == "0/0":
                    status = "FAIL_GT"

                gl_values = gt.get("GL", "").split(",")
                gt_values = gt.get("GT", "").split("/")
                ao_values = (gt.get("AO") or "0").split(",")

                dp = gt.get("DP", "")

                gt0 = int(gt_values[0])
                gt1 = int(gt_values[1])
                idx = int((gt1 * (gt1 + 1) / 2) + gt0)

                depth[sample_type] = dp
                altobs[sample_type] = ao_values
                genotype[sample_type] = gt_values
                likelihood[sample_type] = gl_values
                gl_idx[sample_type] = idx

            lod_norm = leading_float(likelihood["N"][gl_idx["N"]]) - leading_float(likelihood["N"][gl_idx["T"]])
            lod_tumor = leading_float(likelihood["T"][gl_idx["T"]]) - leading_float(likelihood["T"][gl_idx["N"]])
            dqual = lod_tumor + lod_norm

            if genotype["T"][0] == genotype["T"][1] and (
                genotype["T"][1] == genotype["N"][0] or genotype["T"][1] == genotype["N"][1]
            ):
                status = "LOH"

            if dqual < ssc_threshold:
                status = "FAIL_QUAL"

            if lod_norm < lod_threshold or lod_tumor < lod_threshold:
                status = "FAIL_LOD"

            if int(genotype["T"][1]) != int(genotype["N"][0]) and int(genotype["T"][1]) != int(genotype["N"][1]):
                talt = genotype["T"][1]
            else:
                talt = genotype["T"][0]
                if talt == "0":
                    status = "WARN_NOVAR"

            if leading_float(depth["N"]) > 0 and leading_float(depth["T"]) > 0:
                talt_idx = int(talt) - 1
                nvaf = leading_float(altobs["N"][talt_idx]) / leading_float(depth["N"])
                tvaf = leading_float(altobs["T"][talt_idx]) / leading_float(depth["T"])
                if nvaf > 0 and (tvaf / nvaf < min_vaf_ratio):
                    status = "FAIL_NVAF"

            statuses = [status]
            if leading_float(var["INFO"].get("MQM")) <= min_mapping_quality:
                statuses.append("WARN_MQ")

            var["FILTER"] += ";" + ";".join(statuses)
            add_info(var, "SSC", format_compact_number(dqual))

            vcfstr(var, out_fh)


def build_parser():
    """Build the command-line parser.

    Returns:
        Configured argument parser.
    """
    parser = ArgumentParser(description="Filter paired tumor/normal FreeBayes VCF records.")
    parser.add_argument("--vcf", dest="vcf_file", required=True, help="Input FreeBayes VCF file.")
    parser.add_argument("--tumor", dest="tumor_sample", required=True, help="Tumor sample name.")
    parser.add_argument("--normal", dest="normal_sample", required=True, help="Normal sample name.")
    parser.add_argument("--out", dest="out_file", required=True, help="Output filtered VCF file.")
    parser.add_argument(
        "--ssc-threshold", type=float, default=SSC_THRES, help="Minimum somatic score before FAIL_QUAL."
    )
    parser.add_argument(
        "--lod-threshold", type=float, default=LOD_THRES, help="Minimum tumor and normal LOD before FAIL_LOD."
    )
    parser.add_argument("--min-vaf-ratio", type=float, default=MIN_VAF_RATIO, help="Minimum tumor/normal VAF ratio.")
    parser.add_argument(
        "--min-mapping-quality", type=float, default=10, help="MQM value at or below this adds WARN_MQ."
    )
    return parser


def main(argv=None):
    """Run the command-line interface.

    Args:
        argv: Optional command-line argument vector.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args(argv)
    filter_vcf(
        args.vcf_file,
        args.tumor_sample,
        args.normal_sample,
        args.out_file,
        args.ssc_threshold,
        args.lod_threshold,
        args.min_vaf_ratio,
        args.min_mapping_quality,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
