#!/usr/bin/env python3
"""Filter Manta SV records."""

import re
from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import leading_float


def _manta_passes_af(var, sample_id, min_af):
    """Return whether a Manta record has enough PR/SR support."""
    for gt in var["GT"]:
        if gt["_sample_id"] != sample_id:
            continue
        pr = [0.0, 0.0]
        sr = [0.0, 0.0]
        if gt.get("PR"):
            pr = [leading_float(item) for item in gt["PR"].split(",")[:2]]
        if gt.get("SR"):
            sr = [leading_float(item) for item in gt["SR"].split(",")[:2]]
        if (pr[0] == 0 and sr[0] == 0) or ((pr[0] + pr[1]) < 20 and (sr[0] + sr[1]) < 20):
            return False
        af_pr = pr[1] / (pr[0] + pr[1]) if pr[0] + pr[1] else 0
        af_sr = sr[1] / (sr[0] + sr[1]) if sr[0] + sr[1] else 0
        if af_sr < min_af and af_pr < min_af:
            return False
    return True


def filter_manta(vcf_file, sample_id, min_af, out_prefix=None, filtered_out=None, bnd_out=None):
    """Split and filter Manta BND and non-BND records."""
    prefix = out_prefix or sample_id
    bnd_path = bnd_out or f"{prefix}_manta_bnd_filtered.vcf"
    filtered_path = filtered_out or f"{prefix}_manta_filtered.vcf"
    with vcf.VCFReader(vcf_file) as reader, open(bnd_path, "w") as bnd_fh, open(filtered_path, "w") as filtered_fh:
        for line in reader.meta_header_lines:
            if line.startswith("##contig") and not re.search(r"ID=[0-9XYM]{1,2},", line):
                continue
            print(line, file=bnd_fh)
            print(line, file=filtered_fh)
        print(reader.column_header_line, file=bnd_fh)
        print(reader.column_header_line, file=filtered_fh)

        for var in reader:
            if not _manta_passes_af(var, sample_id, leading_float(min_af)):
                continue
            if var["INFO"].get("SVTYPE") == "BND":
                print(var["vcf_str"], file=bnd_fh)
            else:
                print(var["vcf_str"], file=filtered_fh)
    return bnd_path, filtered_path


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input Manta VCF.")
    parser.add_argument("--sample-id", "--id", dest="sample_id", required=True, help="Sample ID.")
    parser.add_argument("--af", required=True, type=float, help="Minimum PR/SR allele fraction.")
    parser.add_argument("--out-prefix", help="Output prefix. Defaults to sample ID.")
    parser.add_argument("--filtered-out", help="Output VCF for non-BND Manta records.")
    parser.add_argument("--bnd-out", help="Output VCF for BND Manta records.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    filter_manta(args.vcf, args.sample_id, args.af, args.out_prefix, args.filtered_out, args.bnd_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
