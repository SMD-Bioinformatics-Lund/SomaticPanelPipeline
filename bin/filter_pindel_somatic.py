#!/usr/bin/env python3
"""Remove Pindel VCF records that are homozygous reference in all samples."""

from argparse import ArgumentParser
from pathlib import Path


def filter_vcf(vcf_file, out_file):
    """Filter Pindel VCF records.

    Args:
        vcf_file: Input VCF path.
        out_file: Output VCF path.
    """
    input_path = Path(vcf_file)
    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise FileNotFoundError(f"{vcf_file} does not exist or is empty")

    with input_path.open() as in_fh, open(out_file, "w") as out_fh:
        for line in in_fh:
            if line.startswith("#"):
                out_fh.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if any(not sample.startswith("0/0") for sample in fields[9:]):
                out_fh.write(line)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Filter Pindel VCF records.")
    parser.add_argument("--vcf", required=True, help="Input Pindel VCF file.")
    parser.add_argument("--out", required=True, help="Output filtered VCF file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    filter_vcf(args.vcf, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
