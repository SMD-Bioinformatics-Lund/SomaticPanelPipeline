#!/usr/bin/env python3
"""Convert GATK segment table to fold-change table."""

from argparse import ArgumentParser

from vcf_pipeline_utils import leading_float, open_output


def gatk_segments(cr_file, purity, out_file):
    """Convert GATK segment table to fold-change table."""
    with open(cr_file) as in_fh, open_output(out_file) as out_fh:
        for line in in_fh:
            if line.startswith("@") or line.startswith("C"):
                continue
            fields = line.rstrip("\n").split("\t")
            foldchange = (2 ** leading_float(fields[4])) / leading_float(purity)
            out_fh.write(f"{fields[0]}\t{fields[4]}\t{foldchange}\t{fields[5]}\n")


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cr", required=True, help="Input copy-ratio segment table.")
    parser.add_argument("--baf", help="Accepted for pipeline compatibility; ignored.")
    parser.add_argument("--purity", required=True, help="Tumor purity.")
    parser.add_argument("--out", required=True, help="Output table.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    gatk_segments(args.cr, args.purity, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
