#!/usr/bin/env python3
"""Filter unpaired FreeBayes VCF records for variant review.

The script reads a single-sample FreeBayes VCF, appends filter tags based on
genotype and mapping quality fields, and writes the filtered VCF to a file.
"""

import re
import sys
from argparse import ArgumentParser
from contextlib import nullcontext

import vcf2 as vcf2

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
    """Write the original VCF header plus unpaired filter metadata.

    Args:
        reader: Open VCF reader with parsed header lines.
        out_fh: Writable output file handle.
    """
    for line in reader.meta_header_lines:
        print(line, file=out_fh)

    print('##FILTER=<ID=FAIL_QUAL,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=WARN_NOVAR,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=WARN_MQ,Description="Record fails the filters">', file=out_fh)
    print('##FILTER=<ID=.,Description="Record fails the filters">', file=out_fh)
    if reader.column_header_line:
        print(reader.column_header_line, file=out_fh)


def filter_vcf(vcf_file, out_file, min_mapping_quality=10):
    """Filter an unpaired FreeBayes VCF.

    Args:
        vcf_file: Path to the input VCF file.
        out_file: Path to the output VCF file.
    """
    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        print_header(reader, out_fh)
        for var in reader:
            genotype = {}
            status = "PASS"

            for gt in var["GT"]:
                sample_type = "T"

                if sample_type == "T" and gt.get("GT") == "0/0":
                    status = "FAIL_GT"

                gt_values = gt.get("GT", "").split("/")
                genotype[sample_type] = gt_values

            talt = genotype["T"][0]
            if talt == "0":
                status = "WARN_NOVAR"

            statuses = [status]
            if leading_float(var["INFO"].get("MQM")) <= min_mapping_quality:
                statuses.append("WARN_MQ")

            var["FILTER"] += ";" + ";".join(statuses)

            vcfstr(var, out_fh)


def build_parser():
    """Build the command-line parser.

    Returns:
        Configured argument parser.
    """
    parser = ArgumentParser(description="Filter unpaired FreeBayes VCF records.")
    parser.add_argument("--vcf", dest="vcf_file", required=True, help="Input FreeBayes VCF file.")
    parser.add_argument("--out", dest="out_file", required=True, help="Output filtered VCF file.")
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
    filter_vcf(args.vcf_file, args.out_file, args.min_mapping_quality)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
