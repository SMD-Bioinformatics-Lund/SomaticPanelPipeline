#!/usr/bin/env python3
"""Normalize multi-allelic VEP gnomAD fields in a VCF."""

import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf2  # noqa: E402
from vcf_pipeline_utils import leading_float, open_output, write_filter_header, write_variant  # noqa: E402


FIELDS_TO_FIX = {"gnomADg", "gnomADg_AF", "gnomADg_AF_popmax", "gnomADg_popmax"}


def csq_field_order(reader):
    """Return the CSQ annotation field order from the VCF header.

    Args:
        reader: Open VCFReader instance.

    Returns:
        List of CSQ field names.
    """
    description = reader.meta.get("INFO", {}).get("CSQ", {}).get("Description", "")
    marker = "Consequence annotations from Ensembl VEP. Format: "
    field_str = description.split(marker, 1)[1] if marker in description else description
    fields = field_str.split("|") if field_str else []
    if fields:
        fields[0] = "Allele"
    return fields


def max_ampersand_value(value):
    """Return the largest numeric value from an ampersand-delimited field.

    Args:
        value: Ampersand-delimited numeric text.

    Returns:
        Largest value as compact string.
    """
    maximum = 0.0
    for item in str(value or "").split("&"):
        number = leading_float(item)
        if number > maximum:
            maximum = number
    return f"{maximum:g}"


def fixed_value(key, value):
    """Return the normalized value for one VEP CSQ field.

    Args:
        key: CSQ field name.
        value: CSQ field value.

    Returns:
        Normalized field value.
    """
    if key not in FIELDS_TO_FIX:
        if key == "Consequence" and isinstance(value, list):
            return "&".join(value)
        return value

    values = str(value or "").split("&")
    if len(values) < 2:
        return value
    if key in ("gnomADg", "gnomADg_popmax"):
        return values[0]
    return max_ampersand_value(value)


def fix_vcf(vcf_file, out_file):
    """Normalize VEP gnomAD annotations in a VCF.

    Args:
        vcf_file: Input VCF path.
        out_file: Output VCF path.
    """
    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        write_filter_header(reader, out_fh, [])
        fields = csq_field_order(reader)
        for var in reader:
            csq_items = []
            for transcript in var["INFO"].get("CSQ", []):
                values = [str(fixed_value(key, transcript.get(key, ""))) for key in fields]
                csq_items.append("|".join(values))
            var["INFO"]["_CSQ_str"] = ",".join(csq_items)
            write_variant(var, out_fh, use_original_csq=True, trailing_sample_tab=False)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Normalize multi-allelic VEP gnomAD fields.")
    parser.add_argument("--vcf", required=True, help="Input VCF file.")
    parser.add_argument("--out", required=True, help="Output VCF file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    fix_vcf(args.vcf, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
