"""Germline-related VCF utilities."""

from __future__ import annotations

import json
from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import add_info, leading_float, open_output, write_filter_header, write_variant


def max_consequence_score(consequences, rank):
    """Return the highest consequence score from an assay definition."""
    scores = rank.get("consequence_score", {})
    return max((scores.get(item, 0) for item in consequences), default=0)


def passes_assay_rank(var, assay):
    """Return whether a germline candidate passes assay inclusion scoring."""
    score = 0
    clinvar_score = 0
    if var["INFO"].get("CLNSIG") and assay.get("clinvar"):
        for value in str(var["INFO"]["CLNSIG"]).split("/"):
            clinvar_score = max(clinvar_score, assay["clinvar"].get(value, 0))
    score += clinvar_score

    max_score = 0
    gnomad = 0
    for tx in var["INFO"].get("CSQ", []):
        if tx.get("Consequence") and assay.get("consequence_cutoff") and assay.get("consequence_score"):
            max_score = max(max_score, max_consequence_score(tx["Consequence"], assay))
        if (tx.get("gnomADg_AF") or tx.get("gnomAD_AF")) and assay.get("gnomad_cutoff"):
            genome = leading_float(str(tx.get("gnomADg_AF", "")).split(",")[0])
            exome = leading_float(str(tx.get("gnomAD_AF", "")).split(",")[0])
            if (genome or exome) < leading_float(assay["gnomad_cutoff"]):
                gnomad = 1

    if assay.get("consequence_cutoff") and max_score >= leading_float(assay["consequence_cutoff"]):
        score += 1
    if gnomad >= 1:
        score += 1
    return score >= leading_float(assay.get("inclusion_score"))


def mark_germlines(
    vcf_file,
    tumor_id,
    normal_id,
    assay_file,
    out_file,
    min_normal_vaf=0.35,
    min_tumor_vaf=0.45,
    min_depth=100,
):
    """Add GERMLINE and GERMLINE_RISK filters to a VCF."""
    with open(assay_file) as fh:
        assay = json.load(fh)
    genes = {gene: True for gene in assay.get("genes", [])}

    headers = [
        '##FILTER=<ID=GERMLINE,Description="Germline variant, detected in normal sample">',
        '##FILTER=<ID=GERMLINE_RISK,Description="Potential germline variant, from tumor sample">',
    ]
    with vcf.VCFReader(vcf_file) as reader, open_output(out_file) as out_fh:
        write_filter_header(reader, out_fh, headers)
        for var in reader:
            in_gene = any(genes.get(tx.get("SYMBOL")) or genes.get("ALL_GENES") for tx in var["INFO"].get("CSQ", []))
            germline = False
            germline_risk = False
            for gt in var["GT"]:
                if normal_id and gt["_sample_id"] == normal_id:
                    germline = (
                        leading_float(gt.get("VAF")) > min_normal_vaf
                        and leading_float(gt.get("DP")) > min_depth
                        and in_gene
                    )
                if tumor_id and gt["_sample_id"] == tumor_id:
                    germline_risk = (
                        leading_float(gt.get("VAF")) > min_tumor_vaf and leading_float(gt.get("DP")) > min_depth
                    )

            if assay.get("inclusion_score") and germline:
                germline = passes_assay_rank(var, assay)

            if germline:
                filters = ["GERMLINE"] + [item for item in var["FILTER"].split(";") if item != "FAIL_NVAF"]
                var["FILTER"] = ";".join(filters)
            elif germline_risk and not normal_id:
                var["FILTER"] = ";".join(var["FILTER"].split(";") + ["GERMLINE_RISK"])
            write_variant(var, out_fh, use_original_csq=True, trailing_sample_tab=False)


def germline_for_cnvkit(vcf_file, out_file, min_gnomad_af=0.05, min_variant_depth=50):
    """Keep common germline variants suitable for CNVkit BAF input."""
    with vcf.VCFReader(vcf_file) as reader, open_output(out_file) as out_fh:
        for line in reader.meta_header_lines:
            if line.startswith("##INFO") and "AF,Number=A" not in line:
                print(line, file=out_fh)
                print(
                    '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">',
                    file=out_fh,
                )
            elif "ID=VAF" not in line:
                print(line, file=out_fh)
        print(reader.column_header_line, file=out_fh)
        for var in reader:
            vaf_value = ""
            vd = 0
            for gt in var["GT"]:
                vaf_value = gt.get("VAF", "")
                vd = leading_float(gt.get("VD"))
            gnomad = 0
            if var["INFO"].get("CSQ") and var["INFO"]["CSQ"][0].get("gnomADg_AF"):
                gnomad = leading_float(var["INFO"]["CSQ"][0].get("gnomADg_AF"))
            if not (gnomad >= min_gnomad_af and vd >= min_variant_depth):
                continue
            add_info(var, "AF", vaf_value)
            write_variant(var, out_fh, use_original_csq=True)


def build_mark_parser():
    """Build the mark-germlines parser."""
    parser = ArgumentParser(description="Mark likely germline variants in a VCF.")
    parser.add_argument("--vcf", required=True, help="Input VCF file.")
    parser.add_argument("--tumor-id", help="Tumor sample ID.")
    parser.add_argument("--normal-id", help="Normal sample ID.")
    parser.add_argument("--assay", required=True, help="Assay JSON file.")
    parser.add_argument("--out", required=True, help="Output VCF file.")
    parser.add_argument("--min-normal-vaf", type=float, default=0.35)
    parser.add_argument("--min-tumor-vaf", type=float, default=0.45)
    parser.add_argument("--min-depth", type=float, default=100)
    return parser


def main_mark(argv=None):
    """Run mark-germlines CLI."""
    args = build_mark_parser().parse_args(argv)
    mark_germlines(
        args.vcf,
        args.tumor_id,
        args.normal_id,
        args.assay,
        args.out,
        args.min_normal_vaf,
        args.min_tumor_vaf,
        args.min_depth,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main_mark())
