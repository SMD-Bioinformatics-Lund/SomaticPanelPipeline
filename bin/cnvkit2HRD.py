#!/usr/bin/env python3
"""Convert CNVkit calls to scarHRD input format."""

from argparse import ArgumentParser

from vcf_pipeline_utils import leading_float, open_output


def cnvkit_to_hrd(input_file, sample_id, ploidy, out_file):
    """Convert CNVkit call table to HRD input format."""
    with open(input_file) as in_fh, open_output(out_file) as out_fh:
        suffix = "\tploidy" if ploidy else ""
        out_fh.write(f"SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn{suffix}\n")
        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith("chromosome") or not line:
                continue
            fields = line.split("\t")
            if fields[0] == "Y" and leading_float(fields[8]) == 0:
                continue
            a_cn = str(leading_float(fields[8]) - 1) if fields[10] == "" or fields[9] == "" else fields[10]
            b_cn = "1" if fields[10] == "" or fields[9] == "" else fields[9]
            row = [sample_id, f"chr{fields[0]}", fields[1], fields[2], fields[8], a_cn, b_cn]
            if ploidy:
                row.append(ploidy)
            out_fh.write("\t".join(row) + "\n")


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input CNVkit call table.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--ploidy", help="Optional ploidy value.")
    parser.add_argument("--out", required=True, help="Output table.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    cnvkit_to_hrd(args.input, args.sample_id, args.ploidy, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
