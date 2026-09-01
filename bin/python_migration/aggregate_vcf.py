#!/usr/bin/env python3
"""Aggregate prioritized caller VCF files into one compact VCF."""

import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf2  # noqa: E402
from vcf_pipeline_utils import add_info, format_compact_number, leading_float, open_output  # noqa: E402


SUPPORTED_CALLERS = ("freebayes", "mutect2", "tnscope", "vardict", "pindel")


def parse_vcf_list(value):
    """Parse and validate a comma-delimited VCF list.

    Args:
        value: Comma-delimited VCF paths.

    Returns:
        List of VCF paths.
    """
    files = []
    for item in value.split(","):
        path = Path(item)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"VCF does not exist: {item}")
        files.append(item)
    return files


def variant_id(var):
    """Return the aggregate key for a variant."""
    return f"{var['CHROM']}_{var['POS']}_{var['REF']}_{var['ALT']}"


def which_variantcaller(meta):
    """Identify the originating variant caller from VCF metadata.

    Args:
        meta: Parsed VCF metadata.

    Returns:
        Caller name.
    """
    source = str(meta.get("source", ""))
    if "freeBayes" in source:
        return "freebayes"
    if "Mutect2" in source:
        return "mutect2"
    if "pindel" in source:
        return "pindel"
    if meta.get("SentieonCommandLine.TNscope"):
        return "tnscope"
    if meta.get("INFO", {}).get("MSILEN"):
        return "vardict"
    if meta.get("INFO", {}).get("ASSESS"):
        return "MELT"
    return "unknown"


def is_weird_freebayes(var):
    """Return whether a freeBayes record lacks AO in the first sample."""
    if not var["GT"]:
        return True
    return "AO" not in var["GT"][0] and "AO" not in var["INFO"]


def summarize_filters(filters):
    """Summarize caller-specific FILTER values for an aggregated variant.

    Args:
        filters: Iterable of FILTER values.

    Returns:
        Aggregated FILTER string.
    """
    non_pass = []
    pass_filter = True
    for item in filters:
        if item not in ("PASS", "."):
            non_pass.append(item)
        if "FAIL" in item:
            pass_filter = False
    if pass_filter:
        non_pass.insert(0, "PASS")
    return ";".join(non_pass) if non_pass else "PASS"


def add_gt(var, sample, key, value):
    """Add a genotype FORMAT value for a sample.

    Args:
        var: Parsed variant dictionary.
        sample: Sample ID.
        key: FORMAT key.
        value: FORMAT value.
    """
    if not any(existing.startswith(key) for existing in var["FORMAT"]):
        var["FORMAT"].append(key)
    for gt in var["GT"]:
        if gt["_sample_id"] == sample:
            gt[key] = str(value)


def split_ad(value):
    """Return reference and alternate depths from an AD field."""
    parts = str(value or "").split(",")
    ref_depth = leading_float(parts[0]) if len(parts) > 0 else 0.0
    alt_depth = leading_float(parts[1]) if len(parts) > 1 else 0.0
    return ref_depth, alt_depth


def split_info_list(value):
    """Split a comma-delimited INFO value into scalar values."""
    if value is None:
        return []
    return [item for item in str(value).split(",") if item != ""]


def freebayes_old_variant_index(var):
    """Return the allele index for VT-decomposed FreeBayes INFO arrays.

    VT keeps the original multiallelic allele strings in ``OLD_MULTIALLELIC``
    and the selected original allele in ``OLD_VARIANT``. FreeBayes stores
    allele-specific depths such as AO as INFO arrays after decomposition, so
    using the first value is wrong for non-first ALT alleles.
    """
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


def select_freebayes_allele_value(var, gt, key):
    """Return the allele-specific FreeBayes FORMAT/INFO value for a decomposed ALT."""
    values = split_info_list(gt.get(key))
    if not values:
        values = split_info_list(var["INFO"].get(key))
    if not values:
        return None
    idx = freebayes_old_variant_index(var)
    if idx < len(values):
        return values[idx]
    return values[0]


def perl_truthy(value):
    """Return truthiness matching Perl scalar checks for simple strings."""
    return value not in (None, "", "0")


def fix_gt(var, caller):
    """Normalize genotype FORMAT fields for aggregate output.

    Args:
        var: Parsed variant dictionary.
        caller: Variant caller name.
    """
    var["FORMAT"] = []
    if caller in ("mutect2", "tnscope", "vardict", "pindel"):
        for gt in var["GT"]:
            ref_depth, alt_depth = split_ad(gt.get("AD"))
            total_depth = ref_depth + alt_depth
            calculated_af = 0.0
            if total_depth > 0:
                calculated_af = alt_depth / total_depth
            af = gt.get("AF") if perl_truthy(gt.get("AF")) else format_compact_number(calculated_af)
            add_gt(var, gt["_sample_id"], "GT", gt.get("GT", ""))
            add_gt(var, gt["_sample_id"], "VAF", af)
            add_gt(var, gt["_sample_id"], "VD", format_compact_number(alt_depth))
            add_gt(var, gt["_sample_id"], "DP", format_compact_number(total_depth))
    elif caller == "freebayes":
        for gt in var["GT"]:
            vaf = 0
            alt_depth = 0
            ao = select_freebayes_allele_value(var, gt, "AO")
            dp = gt.get("DP")
            if not perl_truthy(dp) and perl_truthy(var["INFO"].get("DP")):
                dp = var["INFO"].get("DP")
            if perl_truthy(ao) and ao != ".":
                vaf = f"{leading_float(ao) / leading_float(dp):.4f}"
                alt_depth = ao
            add_gt(var, gt["_sample_id"], "GT", gt.get("GT", ""))
            add_gt(var, gt["_sample_id"], "VAF", vaf)
            add_gt(var, gt["_sample_id"], "VD", alt_depth)
            add_gt(var, gt["_sample_id"], "DP", dp or "")
    elif caller == "MELT":
        for gt in var["GT"]:
            add_gt(var, gt["_sample_id"], "GT", gt.get("GT", ""))
            add_gt(var, gt["_sample_id"], "VAF", gt.get("VAF", ""))
            add_gt(var, gt["_sample_id"], "VD", gt.get("VD", ""))
            add_gt(var, gt["_sample_id"], "DP", gt.get("DP", ""))


def aggregate_vcfs(vcf_files):
    """Aggregate variants from prioritized input VCFs.

    Args:
        vcf_files: Input VCF paths in priority order.

    Returns:
        Tuple of aggregated variant dictionary, all observed filters, and first
        reader column header line.
    """
    aggregated = {}
    all_filters = {}
    filters_by_variant = {}
    first_column_header = ""

    for vcf_file in vcf_files:
        with vcf2.VCFReader(file=vcf_file) as reader:
            if not first_column_header:
                first_column_header = reader.column_header_line
            caller = which_variantcaller(reader.meta)
            for var in reader:
                simple_id = variant_id(var)
                melt_info_keep = var["INFO"].get("SCOUT_CUSTOM", 0)

                if var.get("FILTER"):
                    for item in var["FILTER"].split(";"):
                        filters_by_variant.setdefault(simple_id, {})[item] = True
                        all_filters[item] = True

                if caller == "freebayes" and is_weird_freebayes(var):
                    continue

                if simple_id in aggregated:
                    aggregated[simple_id]["INFO"]["variant_callers"] += f"|{caller}"
                else:
                    fix_gt(var, caller)
                    var["INFO"] = {}
                    var["INFO_order"] = []
                    add_info(var, "variant_callers", caller)
                    if caller == "MELT":
                        add_info(var, "custom", melt_info_keep)
                    aggregated[simple_id] = var

    for simple_id, var in aggregated.items():
        var["FILTER"] = summarize_filters(filters_by_variant.get(simple_id, {}).keys())
    return aggregated, list(all_filters.keys()), first_column_header


def write_header(out_fh, filters, vcf_files, sample_order, first_column_header):
    """Write the aggregate VCF header."""
    print("##fileformat=VCFv4.2", file=out_fh)
    print("##origin=" + ",".join(vcf_files), file=out_fh)
    print('##INFO=<ID=variant_callers,Number=.,Type=String,Description="List of variant callers which detected the variant">', file=out_fh)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=out_fh)
    print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', file=out_fh)
    print('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="ALT allele observation fraction">', file=out_fh)
    print('##FORMAT=<ID=VD,Number=1,Type=Integer,Description="ALT allele observation count">', file=out_fh)
    for item in filters:
        print(f'##FILTER=<ID={item},Description="{item}">', file=out_fh)
    if sample_order:
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_order), file=out_fh)
    elif first_column_header:
        print(first_column_header, file=out_fh)


def serialize_aggregate_variant(var, sample_order):
    """Serialize an aggregated variant.

    Args:
        var: Parsed variant dictionary.
        sample_order: Optional sample output order.

    Returns:
        VCF record line.
    """
    info = ";".join(f"{key}={var['INFO'][key]}" for key in var["INFO_order"])
    fields = [
        var["CHROM"],
        var["POS"],
        var["ID"],
        var["REF"],
        var["ALT"],
        var["QUAL"],
        var["FILTER"],
        info,
        ":".join(var["FORMAT"]),
    ]
    if sample_order:
        order = {sample: idx for idx, sample in enumerate(sample_order)}
    else:
        order = {gt["_sample_id"]: idx for idx, gt in enumerate(var["GT"])}
    for gt in sorted(var["GT"], key=lambda item: order[item["_sample_id"]]):
        fields.append(":".join(gt.get(key, "") for key in var["FORMAT"]))
    return "\t".join(fields)


def aggregate(vcf_arg, sample_order_arg, out_file):
    """Aggregate VCFs and write output.

    Args:
        vcf_arg: Comma-delimited input VCF paths.
        sample_order_arg: Optional comma-delimited sample order.
        out_file: Output VCF path.
    """
    vcf_files = parse_vcf_list(vcf_arg)
    sample_order = sample_order_arg.split(",") if sample_order_arg else []
    aggregated, filters, first_column_header = aggregate_vcfs(vcf_files)
    with open_output(out_file) as out_fh:
        write_header(out_fh, filters, vcf_files, sample_order, first_column_header)
        for var in aggregated.values():
            out_fh.write(serialize_aggregate_variant(var, sample_order) + "\n")


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Aggregate prioritized caller VCF files.")
    parser.add_argument("--vcf", required=True, help="Comma-delimited VCF files in priority order.")
    parser.add_argument("--sample-order", help="Comma-delimited output sample order.")
    parser.add_argument("--out", required=True, help="Output aggregate VCF file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    aggregate(args.vcf, args.sample_order, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
