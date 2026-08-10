#!/usr/bin/env python3
"""Merge adjacent tumor GATK/CNV SV records."""

from mergeGATK import merge_gatk
from argparse import ArgumentParser


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--out", required=True, help="Output VCF.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    merge_gatk(args.vcf, args.out, tumor=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
