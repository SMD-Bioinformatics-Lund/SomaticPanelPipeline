#!/usr/bin/env python3
"""Filter CNVkit TSV calls."""

from argparse import ArgumentParser

from vcf_pipeline_utils import leading_float, open_output


def filter_cnvkit(tsv_file, coverage, distance, out_file):
    """Filter CNVkit TSV calls."""
    cov = leading_float(coverage)
    max_dist = leading_float(distance)
    with open(tsv_file) as in_fh, open_output(out_file) as out_fh:
        for idx, line in enumerate(in_fh, start=1):
            fields = line.split("\t")
            if idx == 1:
                out_fh.write(line)
                continue
            if fields[3] == "-":
                continue
            if leading_float(fields[6]) == 2:
                continue
            if leading_float(fields[2]) - leading_float(fields[1]) > max_dist:
                continue
            if leading_float(fields[6]) == 1 and cov and leading_float(fields[9]) / cov <= 0.10:
                continue
            out_fh.write("\t".join(fields))


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True, help="Input CNVkit TSV.")
    parser.add_argument("--coverage", required=True, help="Coverage value.")
    parser.add_argument("--distance", required=True, help="Maximum segment length.")
    parser.add_argument("--out", required=True, help="Output TSV.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    filter_cnvkit(args.tsv, args.coverage, args.distance, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
