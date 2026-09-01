#!/usr/bin/env python3
"""Germline-related VCF utilities."""

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
        af_header_written = False
        for line in reader.meta_header_lines:
            if line.startswith("##INFO") and "AF,Number=A" in line:
                af_header_written = True
                print(line, file=out_fh)
            elif line.startswith("##INFO"):
                print(line, file=out_fh)
                if not af_header_written:
                    print(
                        '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">',
                        file=out_fh,
                    )
                    af_header_written = True
            elif "ID=VAF" not in line:
                print(line, file=out_fh)
        print(reader.column_header_line, file=out_fh)
        csq_fields = csq_field_order(reader)
        for var in reader:
            vaf_value = ""
            vd = 0
            for gt in var["GT"]:
                vaf_value = gt.get("VAF", "")
                vd = leading_float(gt.get("VD"))
            gnomad = first_csq_float(var, "gnomADg_AF", csq_fields)
            if not (gnomad >= min_gnomad_af and vd >= min_variant_depth):
                continue
            add_info(var, "AF", vaf_value)
            write_variant(var, out_fh, use_original_csq=True)


def first_csq_float(var, field_name, csq_fields):
    """Return a numeric CSQ field from the first transcript.

    The Perl helper reads ``INFO/CSQ[0]`` after vcf2 has parsed the VEP CSQ
    string. Keep the same behavior, but fall back to the original CSQ string
    when the lightweight Python parser cannot expose the field cleanly.
    """
    csq = var["INFO"].get("CSQ") or []
    if csq and isinstance(csq[0], dict):
        value = csq[0].get(field_name, "")
        if isinstance(value, list):
            value = value[0] if value else ""
        parsed = max_ampersand_float(value)
        if parsed:
            return parsed

    csq_text = var["INFO"].get("_CSQ_str", "")
    if not csq_text or field_name not in csq_fields:
        return 0.0
    field_idx = csq_fields.index(field_name)
    first_transcript = csq_text.split(",", 1)[0].split("|")
    if field_idx >= len(first_transcript):
        return 0.0
    return max_ampersand_float(first_transcript[field_idx])


def csq_field_order(reader):
    """Return CSQ field order from the VCF header."""
    description = reader.meta.get("INFO", {}).get("CSQ", {}).get("Description", "")
    marker = "Consequence annotations from Ensembl VEP. Format: "
    if marker in description:
        description = description.split(marker, 1)[1]
    return description.split("|") if description else []


def max_ampersand_float(value):
    """Return the largest leading float in an ampersand-delimited CSQ value."""
    best = 0.0
    for item in str(value or "").split("&"):
        best = max(best, leading_float(item))
    return best


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


def build_parser():
    parser = ArgumentParser(description="Create a CNVkit-compatible germline VCF.")
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--out", required=True, help="Output VCF.")
    parser.add_argument("--min-gnomad-af", type=float, default=0.05)
    parser.add_argument("--min-variant-depth", type=float, default=50)
    return parser

def main(argv=None):
    args = build_parser().parse_args(argv)
    germline_for_cnvkit(args.vcf, args.out, args.min_gnomad_af, args.min_variant_depth)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
