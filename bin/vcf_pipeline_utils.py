#!/usr/bin/env python3
"""Shared helpers for SomaticPanelPipeline VCF command-line filters."""

import re
import sys
from contextlib import nullcontext


LEADING_FLOAT_RE = re.compile(
    r"^[\s]*([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)"
)


def leading_float(value):
    """Convert a scalar-like value to a number using leading numeric text.

    Args:
        value: Input scalar value.

    Returns:
        Floating-point numeric value. Empty or non-numeric values become 0.0.
    """
    if value is None or value == "":
        return 0.0
    match = LEADING_FLOAT_RE.match(str(value))
    if match:
        return float(match.group(1))
    return 0.0


def format_compact_number(value):
    """Format a numeric value without unnecessary trailing zeros.

    Args:
        value: Numeric value to format.

    Returns:
        Compact decimal string.
    """
    return "{:.15g}".format(value)


def open_output(path):
    """Open an output path or stdout.

    Args:
        path: Output file path, or ``-`` for stdout.

    Returns:
        Context manager yielding a writable text handle.
    """
    if path == "-":
        return nullcontext(sys.stdout)
    return open(path, "w")


def add_info(var, key, value):
    """Append an INFO field to a parsed variant.

    Args:
        var: Parsed variant dictionary.
        key: INFO key to append.
        value: INFO value.
    """
    var["INFO_order"].append(key)
    var["INFO"][key] = value


def append_filter_status(var, statuses):
    """Append filtering statuses using the pipeline's VCF FILTER semantics.

    Args:
        var: Parsed variant dictionary.
        statuses: List of FILTER values to add.
    """
    if not statuses:
        return
    if var["FILTER"] in ("PASS", "."):
        var["FILTER"] = ";".join(statuses)
    else:
        var["FILTER"] += ";" + ";".join(statuses)


def write_filter_header(reader, out_fh, extra_header_lines):
    """Write original VCF metadata, extra headers, and the column header.

    Args:
        reader: Open VCFReader instance.
        out_fh: Writable output handle.
        extra_header_lines: Iterable of header lines to insert before #CHROM.
    """
    for line in reader.meta_header_lines:
        print(line, file=out_fh)
    for line in extra_header_lines:
        print(line, file=out_fh)
    if reader.column_header_line:
        print(reader.column_header_line, file=out_fh)


def serialize_variant(var, use_original_csq=False, trailing_sample_tab=True):
    """Serialize a parsed variant dictionary to one VCF line.

    Args:
        var: Parsed variant dictionary.
        use_original_csq: If true, serialize CSQ from ``_CSQ_str`` when present.
        trailing_sample_tab: Whether to keep the historical trailing tab after
            sample fields.

    Returns:
        Serialized VCF line without the final newline unless requested by the
        caller.
    """
    info_fields = []
    for key in var["INFO_order"]:
        value_key = "_CSQ_str" if use_original_csq and key == "CSQ" else key
        info_fields.append(f"{key}={var['INFO'][value_key]}")

    fields = [
        var["CHROM"],
        var["POS"],
        var["ID"],
        var["REF"],
        var["ALT"],
        var["QUAL"],
        var["FILTER"],
        ";".join(info_fields),
        ":".join(var["FORMAT"]),
    ]
    fields.extend(":".join(gt.get(key, "") for key in var["FORMAT"]) for gt in var["GT"])
    line = "\t".join(fields)
    if trailing_sample_tab:
        line += "\t"
    return line


def write_variant(var, out_fh, use_original_csq=False, trailing_sample_tab=True):
    """Write a parsed variant dictionary to a VCF file handle.

    Args:
        var: Parsed variant dictionary.
        out_fh: Writable output handle.
        use_original_csq: If true, serialize CSQ from ``_CSQ_str`` when present.
        trailing_sample_tab: Whether to retain a trailing sample-field tab.
    """
    out_fh.write(
        serialize_variant(
            var,
            use_original_csq=use_original_csq,
            trailing_sample_tab=trailing_sample_tab,
        )
        + "\n"
    )
