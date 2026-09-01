#!/usr/bin/env python3
"""Generate CNVkit coverage and BAF BED tracks."""

from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import format_compact_number, leading_float


def split_info_list(value):
    """Split a comma-delimited VCF scalar into values."""
    if value is None:
        return []
    return [item for item in str(value).split(",") if item != ""]


def old_variant_index(var):
    """Return the VT-decomposed allele index in original multiallelic arrays."""
    old_variant = var["INFO"].get("OLD_VARIANT")
    old_multiallelic = var["INFO"].get("OLD_MULTIALLELIC")
    if not old_variant or not old_multiallelic:
        return 0

    selected_parts = str(old_variant).split(":", 2)
    multi_parts = str(old_multiallelic).split(":", 2)
    if len(selected_parts) != 3 or len(multi_parts) != 3:
        return 0

    selected_alleles = selected_parts[2].split("/")
    multi_alleles = multi_parts[2].split("/")
    if len(selected_alleles) < 2 or len(multi_alleles) < 2:
        return 0

    selected_alt = selected_alleles[-1]
    for idx, allele in enumerate(multi_alleles[1:]):
        if allele == selected_alt:
            return idx
    return 0


def select_allele_value(var, gt, key):
    """Return the allele-specific value from FORMAT first, then INFO."""
    values = split_info_list(gt.get(key))
    if not values:
        values = split_info_list(var["INFO"].get(key))
    if not values:
        return None
    idx = old_variant_index(var)
    if idx < len(values):
        return values[idx]
    return values[0]


def generate_gens_data(cnr_file, vcf_file, sample_id, out_prefix=None, baf_out=None, cov_out=None):
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
                depth_text = gt.get("DP") or "0"
                depth = leading_float(depth_text)
                allele_depth = select_allele_value(var, gt, "AO")
                if depth < 100 or not allele_depth:
                    break
                vaf_value = leading_float(allele_depth) / depth
                baf_data.append(
                    f"{var['CHROM']}\t{int(var['POS']) - 1}\t{var['POS']}\t"
                    f"{format_compact_number(vaf_value)}"
                )

    prefix = out_prefix or sample_id
    cov_bed = cov_out or f"{prefix}.cov.bed"
    baf_bed = baf_out or f"{prefix}.baf.bed"
    for path, data in ((cov_bed, log2_data), (baf_bed, baf_data)):
        with open(path, "w") as out_fh:
            for resolution in "oabcd":
                for row in data:
                    out_fh.write(f"{resolution}_{row}\n")


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cnr", required=True, help="Input CNVkit CNR file.")
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--out-prefix", help="Output prefix. Defaults to sample ID.")
    parser.add_argument("--baf-out", help="Output BAF BED file.")
    parser.add_argument("--cov-out", help="Output coverage BED file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    generate_gens_data(args.cnr, args.vcf, args.sample_id, args.out_prefix, args.baf_out, args.cov_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
