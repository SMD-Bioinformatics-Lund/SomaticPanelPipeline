#!/usr/bin/env python3
"""Merge adjacent GATK SV records."""

import re
from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import leading_float, open_output


def _parse_info_column(info_text):
    """Parse a raw INFO column into a dictionary."""
    info = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "defined"
    return info


def _natural_key(value):
    """Return a sort key close to GNU sort -V for chromosome strings."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def _format_merged_gatk_record(agg, info, tumor=False):
    """Serialize one merged GATK/CNV segment record."""
    fields = agg[:7]
    if tumor:
        info_items = []
        for key in ("END", "SVLEN", "SVTYPE", "FOLD_CHANGE_LOG", "FOLD_CHANGE", "PROBES"):
            values = str(info[key]).split(",")
            if key in {"FOLD_CHANGE_LOG", "FOLD_CHANGE"}:
                value = sum(leading_float(item) for item in values) / len(values)
            elif key == "PROBES":
                value = int(sum(leading_float(item) for item in values))
            else:
                value = info[key]
            info_items.append(f"{key}={value}")
    else:
        info_items = [f"{key}={info[key]}" for key in ("END", "SVLEN", "SVTYPE", "gatkCN")]
    return "\t".join(fields + [";".join(info_items)] + agg[8:])


def _merge_gatk_type(vcf_file, svtype, tumor=False):
    """Merge adjacent GATK/CNV records for one SV type."""
    headers = []
    variants = []
    agg = []
    agg_info = {}
    merge_key = "FOLD_CHANGE" if tumor else "gatkCN"
    with vcf.open_vcf(vcf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                if not tumor and re.search(r"<ID=gatkCN,", line):
                    line = line.replace("Number=1,", "Number=.,")
                headers.append(line)
                continue
            data = line.rstrip("\n").split("\t")
            info = _parse_info_column(data[7])
            if info.get("SVTYPE") != svtype:
                continue
            if agg:
                dist = abs(int(data[1]) - int(leading_float(agg_info["END"])))
                agg_first = str(agg_info[merge_key]).split(",")[0]
                if (
                    data[0] == agg[0]
                    and dist / abs(leading_float(info["SVLEN"])) < 0.1
                    and dist / abs(leading_float(agg_info["SVLEN"])) < 0.1
                    and dist < 50000
                    and abs(leading_float(agg_first) - leading_float(info[merge_key])) < 0.2
                    and info["SVTYPE"] == agg_info["SVTYPE"]
                ):
                    agg[2] += "," + data[2]
                    agg_info["END"] = info["END"]
                    agg_info["SVLEN"] = leading_float(agg_info["SVLEN"]) + leading_float(info["SVLEN"])
                    if tumor:
                        for key in ("FOLD_CHANGE_LOG", "FOLD_CHANGE", "PROBES"):
                            agg_info[key] = f"{agg_info[key]},{info[key]}"
                    else:
                        agg_info["gatkCN"] = f"{agg_info['gatkCN']},{info['gatkCN']}"
                    continue
                variants.append(_format_merged_gatk_record(agg, agg_info, tumor=tumor))
            agg = data
            agg_info = info
    if agg:
        variants.append(_format_merged_gatk_record(agg, agg_info, tumor=tumor))
    return headers, variants


def merge_gatk(vcf_file, out_file, tumor=False):
    """Merge adjacent GATK segment records."""
    headers, deletions = _merge_gatk_type(vcf_file, "DEL", tumor=tumor)
    _headers, duplications = _merge_gatk_type(vcf_file, "DUP", tumor=tumor)
    records = sorted(
        deletions + duplications,
        key=lambda row: (_natural_key(row.split("\t", 1)[0]), int(row.split("\t", 2)[1])),
    )
    with open_output(out_file) as out_fh:
        out_fh.writelines(headers)
        for row in records:
            out_fh.write(row if row.endswith("\n") else row + "\n")


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--out", required=True, help="Output VCF.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    merge_gatk(args.vcf, args.out, tumor=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
