#!/usr/bin/env python3
"""Utilities for reading VCF headers, records, INFO fields, and genotypes.

This module provides a lightweight parser tailored to the VCF structure used by
the SomaticPanelPipeline helper scripts. It preserves field order for INFO and
sample FORMAT output so filtered VCF records can be reconstructed exactly.
"""

import gzip
from collections import defaultdict


def split_without_trailing_empty_fields(value, sep):
    """Split text while dropping trailing empty fields.

    Args:
        value: String to split.
        sep: Literal separator string.

    Returns:
        A list of fields with trailing empty values removed.
    """
    parts = value.split(sep)
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def open_vcf(fn):
    """Open a plain-text or gzip-compressed VCF file for text reading.

    Args:
        fn: Path to a VCF file.

    Returns:
        A text file handle positioned at the start of the file.
    """
    with open(fn, "rb") as fh:
        is_compressed = fh.read(2) == b"\x1f\x8b"

    if is_compressed:
        return gzip.open(fn, "rt")
    return open(fn, "rt")


class VCFReader:
    """Streaming VCF reader with parsed metadata and sample names.

    Attributes:
        file: Input VCF path supplied during initialization.
        fh: Open file handle for streaming variant records.
        meta: Parsed VCF metadata grouped by metadata type.
        head: Column names from the #CHROM header line.
        samples: Sample names from columns 10 onward.
        meta_header_lines: Original metadata header lines beginning with ##.
        column_header_line: Original #CHROM header line.
        header_str: Original header lines joined as one string.
    """

    def __init__(self, file):
        """Initialize a VCF reader and parse its header.

        Args:
            file: Path to a plain-text or gzip-compressed VCF file.
        """
        self.file = file
        self.fh = None
        self.meta = defaultdict(dict)
        self.head = []
        self.samples = []
        self.meta_header_lines = []
        self.column_header_line = ""
        self.header_str = ""
        self._open_and_parse_header()

    def _open_and_parse_header(self):
        """Open the configured VCF file and parse header metadata.

        Returns:
            The initialized reader instance.
        """
        fn = self.file
        self.fh = open_vcf(fn)

        head = []
        full_header = []
        meta_header_lines = []
        column_header_line = ""

        for raw_line in self.fh:
            line = raw_line.rstrip("\n")

            if line.strip() == "":
                continue

            if line.startswith("#"):
                full_header.append(line + "\n")

            if line.startswith("##"):
                meta_header_lines.append(line)
                meta_type, meta, val = parse_metainfo(line)
                if meta_type == "NONE":
                    self.meta[meta] = val
                elif meta_type is not None:
                    self.meta[meta_type][meta["ID"]] = meta
            elif line.startswith("#"):
                column_header_line = line
                line = line[1:]
                head = line.split("\t")
                break

        self.head = head
        self.samples = head[9:]
        self.meta_header_lines = meta_header_lines
        self.column_header_line = column_header_line
        self.header_str = "".join(full_header)
        return self

    def __iter__(self):
        """Return the reader as its own variant iterator.

        Returns:
            The current reader instance.
        """
        return self

    def __enter__(self):
        """Enter a context manager for the open VCF reader.

        Returns:
            The current reader instance.
        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close the VCF file when leaving a context manager."""
        self.close()

    def __next__(self):
        """Read and parse the next variant record.

        Returns:
            A parsed variant dictionary.

        Raises:
            StopIteration: When no more variant records are available.
        """
        variant = self.next_variant()
        if variant is None:
            raise StopIteration
        return variant

    def next_variant(self):
        """Read and parse the next variant record.

        Returns:
            A parsed variant dictionary, or None when no more records are
            available.
        """
        if self.fh is None:
            return None

        line = self.fh.readline()
        if line:
            line = line.rstrip("\n")
            variant = parse_variant(line, self.head, self.meta)
            if variant.get("CHROM"):
                return variant
        return None

    def close(self):
        """Close the underlying VCF file handle."""
        if self.fh is not None:
            self.fh.close()
            self.fh = None


def parse_metainfo(comment):
    """Parse a VCF metadata header line.

    Args:
        comment: A metadata line beginning with ``##``.

    Returns:
        A tuple ``(meta_type, meta, value)``. Structured FORMAT, INFO, SAMPLE,
        and FILTER entries return ``meta_type`` and a metadata dictionary.
        Other key/value metadata returns ``("NONE", key, value)``.
    """
    if comment.startswith("##"):
        comment = comment[2:]
    if "=" in comment:
        meta_type, data = comment.split("=", 1)
    else:
        meta_type, data = comment, ""

    if meta_type in ("FORMAT", "INFO", "SAMPLE", "FILTER"):
        if data.startswith("<") and data.endswith(">"):
            data = data[1:-1]
        pairs, _order = parse_key_value_pairs(data, "=", ",")
        return meta_type, pairs, 0
    if meta_type and data:
        return "NONE", meta_type, data

    return None, None, None


def parse_variant(var_str, head, meta):
    """Parse one VCF variant record.

    Args:
        var_str: Raw tab-delimited VCF record.
        head: Column names from the VCF header.
        meta: Parsed VCF metadata used for optional CSQ decoding.

    Returns:
        A dictionary containing the fixed VCF columns, parsed INFO data,
        FORMAT keys, genotype dictionaries, and the original record string.
    """
    var_data = var_str.split("\t")
    var = {"vcf_str": var_str}

    for idx in range(7):
        var[head[idx]] = var_data[idx] if idx < len(var_data) else ""

    var["INFO"], var["INFO_order"] = parse_info(var_data[7] if len(var_data) > 7 else "")

    if var["INFO"].get("CSQ"):
        var["INFO"]["_CSQ_str"] = var["INFO"]["CSQ"]
        var["INFO"]["CSQ"] = parse_vep_csq(
            var["INFO"]["CSQ"], meta.get("INFO", {}).get("CSQ", {})
        )

    fmt = var_data[8] if len(var_data) > 8 else ""
    var["FORMAT"] = fmt.split(":")

    var["GT"] = []
    for idx in range(9, len(var_data)):
        sample_id = head[idx] if idx < len(head) else ""
        var["GT"].append(parse_sample_genotype(fmt, var_data[idx], sample_id))

    return var


def parse_sample_genotype(format_str, data, sample_id):
    """Parse a sample genotype field using the FORMAT definition.

    Args:
        format_str: Colon-delimited FORMAT string.
        data: Colon-delimited sample data string.
        sample_id: Sample name from the VCF header.

    Returns:
        A dictionary mapping FORMAT keys to sample values, plus ``_sample_id``.
    """
    fmt = format_str.split(":")
    values = data.split(":")

    gt = {key: values[idx] if idx < len(values) else "" for idx, key in enumerate(fmt)}
    gt["_sample_id"] = sample_id
    return gt


def parse_info(info_str):
    """Parse the INFO column from a VCF record.

    Args:
        info_str: Semicolon-delimited INFO string.

    Returns:
        A tuple containing the INFO key/value dictionary and key order.
    """
    return parse_key_value_pairs(info_str, "=", ";")


def parse_vep_csq(csq_value, csq_meta):
    """Parse Ensembl VEP CSQ annotations into transcript dictionaries.

    Args:
        csq_value: Comma-delimited CSQ annotation value from INFO.
        csq_meta: Parsed INFO metadata for the CSQ field.

    Returns:
        A list of transcript annotation dictionaries.
    """
    description = csq_meta.get("Description", "")
    marker = "Consequence annotations from Ensembl VEP. Format: "
    field_str = ""
    if marker in description:
        field_str = description.split(marker, 1)[1]

    field_names = field_str.split("|")
    transcripts = csq_value.split(",")

    data_transcripts = []
    for transcript_CSQ in transcripts:
        values = transcript_CSQ.split("|")
        data = {}
        for idx, field_name in enumerate(field_names):
            value = values[idx] if idx < len(values) else ""
            if field_name == "Consequence":
                data[field_name] = value.split("&")
            else:
                data[field_name] = "" if value in ("", "0") else value
        data_transcripts.append(data)

    return data_transcripts


def parse_key_value_pairs(value, keyval_sep, pair_sep):
    """Parse delimited key/value pairs while preserving key order.

    Args:
        value: Input string containing delimited pairs.
        keyval_sep: Separator between each key and value.
        pair_sep: Separator between pairs.

    Returns:
        A tuple containing the parsed dictionary and ordered keys. Entries
        without a key/value separator are treated as boolean flags with value 1.
    """
    pair_str = split_without_trailing_empty_fields(value, pair_sep)
    pairs = {}
    order = []

    for pair in pair_str:
        if keyval_sep in pair:
            parts = pair.split(keyval_sep)
            key = parts[0]
            val = parts[1] if len(parts) > 1 else ""
            pairs[key] = val
            order.append(key)
        else:
            pairs[pair] = 1

    return pairs, order


__all__ = [
    "VCFReader",
    "parse_metainfo",
    "parse_variant",
    "parse_sample_genotype",
    "parse_info",
]
