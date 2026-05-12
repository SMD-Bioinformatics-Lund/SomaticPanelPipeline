#!/usr/bin/env python3
"""Generate compressed CNVkit coverage and BAF BED tracks."""

import subprocess
from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import leading_float


def generate_gens_data(cnr_file, vcf_file, sample_id, out_prefix):
    """Generate CNVkit coverage and BAF tracks.

    Args:
        cnr_file: Input CNVkit ``.cnr`` file.
        vcf_file: Input VCF with sample genotypes.
        sample_id: Sample name used to select the VCF genotype column.
        out_prefix: Output prefix. Defaults to ``sample_id`` when empty.
    """
    log2_data = []
    with open(cnr_file) as cnr:
        next(cnr, None)
        for line in cnr:
            chrom, start, end, gene, _depth, log2_value, *_rest = line.rstrip(
                "\n"
            ).split("\t")
            if gene == "Antitarget":
                continue
            midpoint = int(start) + int((int(end) - int(start)) / 2)
            log2_data.append(f"{chrom}\t{midpoint - 1}\t{midpoint}\t{log2_value}")

    baf_data = []
    with vcf.VCFReader(vcf_file) as reader:
        for var in reader:
            for gt in var["GT"]:
                if gt["_sample_id"] != sample_id:
                    continue
                depth = leading_float(gt.get("DP"))
                if depth < 100 or not gt.get("AO"):
                    break
                vaf_value = leading_float(gt.get("AO")) / depth
                baf_data.append(
                    f"{var['CHROM']}\t{int(var['POS']) - 1}\t{var['POS']}\t{vaf_value}"
                )

    prefix = out_prefix or sample_id
    cov_bed = f"{prefix}.cov.bed"
    baf_bed = f"{prefix}.baf.bed"
    for path, data in ((cov_bed, log2_data), (baf_bed, baf_data)):
        with open(path, "w") as out_fh:
            for resolution in "oabcd":
                for row in data:
                    out_fh.write(f"{resolution}_{row}\n")
        subprocess.run(["bgzip", "-f", path], check=True)
        subprocess.run(["tabix", f"{path}.gz"], check=True)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cnr", required=True, help="Input CNVkit CNR file.")
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--out-prefix", help="Output prefix. Defaults to sample ID.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    generate_gens_data(args.cnr, args.vcf, args.sample_id, args.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
