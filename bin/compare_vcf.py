#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter
from pathlib import Path


def open_vcf(path):
    """Open plain VCF or gzipped VCF."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info_flags(info_field):
    """
    Extract INFO flags.

    Example:
      DP=100;SOMATIC;AF=0.23

    SOMATIC is treated as a flag because it has no '='.
    """
    if info_field in {".", ""}:
        return set()

    return {item for item in info_field.split(";") if item and "=" not in item}


def parse_vcf(path):
    """
    Parse minimal VCF fields.

    Variant identity:
      CHROM, POS, REF, ALT

    Multi-allelic ALT entries are split into separate variants.
    """
    variants = {}

    stats = {
        "file": str(path),
        "records": 0,
        "variants": 0,
        "pass_records": 0,
        "filtered_records": 0,
        "filter_counts": Counter(),
        "info_flag_counts": Counter(),
        "chrom_counts": Counter(),
    }

    with open_vcf(path) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 8:
                raise ValueError(f"Malformed VCF line in {path}:\n{line}")

            chrom, pos, _id, ref, alt_field, qual, filt, info = fields[:8]

            stats["records"] += 1
            stats["chrom_counts"][chrom] += 1

            if filt in {"PASS", "."}:
                stats["pass_records"] += 1
            else:
                stats["filtered_records"] += 1

            for f in filt.split(";"):
                if f:
                    stats["filter_counts"][f] += 1

            info_flags = parse_info_flags(info)

            for flag in info_flags:
                stats["info_flag_counts"][flag] += 1

            for alt in alt_field.split(","):
                key = (chrom, int(pos), ref, alt)

                variants[key] = {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "filter": filt,
                    "info_flags": info_flags,
                }

                stats["variants"] += 1

    return variants, stats


def write_variant_table(path, variants):
    """Write variant table to TSV."""
    with open(path, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tFILTER\tINFO_FLAGS\n")

        for key in sorted(variants):
            v = variants[key]
            flags = ",".join(sorted(v["info_flags"])) if v["info_flags"] else "."

            out.write(
                f"{v['chrom']}\t"
                f"{v['pos']}\t"
                f"{v['ref']}\t"
                f"{v['alt']}\t"
                f"{v['filter']}\t"
                f"{flags}\n"
            )


def write_stats(path, stats1, stats2, shared_count, unique1_count, unique2_count):
    """Write basic comparison stats."""
    with open(path, "w") as out:
        out.write("Metric\tVCF1\tVCF2\n")
        out.write(f"File\t{stats1['file']}\t{stats2['file']}\n")
        out.write(f"Records\t{stats1['records']}\t{stats2['records']}\n")
        out.write(f"Variants\t{stats1['variants']}\t{stats2['variants']}\n")
        out.write(f"PASS_records\t{stats1['pass_records']}\t{stats2['pass_records']}\n")
        out.write(
            f"Filtered_records\t{stats1['filtered_records']}\t{stats2['filtered_records']}\n"
        )
        out.write(f"Shared_variants\t{shared_count}\t{shared_count}\n")
        out.write(f"Unique_variants\t{unique1_count}\t{unique2_count}\n")

        out.write("\nFILTER_COUNTS\n")
        out.write("FILTER\tVCF1\tVCF2\n")

        all_filters = set(stats1["filter_counts"]) | set(stats2["filter_counts"])

        for filt in sorted(all_filters):
            out.write(
                f"{filt}\t"
                f"{stats1['filter_counts'].get(filt, 0)}\t"
                f"{stats2['filter_counts'].get(filt, 0)}\n"
            )

        out.write("\nINFO_FLAG_COUNTS\n")
        out.write("FLAG\tVCF1\tVCF2\n")

        all_flags = set(stats1["info_flag_counts"]) | set(stats2["info_flag_counts"])

        for flag in sorted(all_flags):
            out.write(
                f"{flag}\t"
                f"{stats1['info_flag_counts'].get(flag, 0)}\t"
                f"{stats2['info_flag_counts'].get(flag, 0)}\n"
            )

        out.write("\nCHROM_COUNTS\n")
        out.write("CHROM\tVCF1\tVCF2\n")

        all_chroms = set(stats1["chrom_counts"]) | set(stats2["chrom_counts"])

        for chrom in sorted(all_chroms):
            out.write(
                f"{chrom}\t"
                f"{stats1['chrom_counts'].get(chrom, 0)}\t"
                f"{stats2['chrom_counts'].get(chrom, 0)}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Compare two VCF files and report only non-shared variants."
    )

    parser.add_argument("vcf1", help="First VCF file. Supports .vcf and .vcf.gz")
    parser.add_argument("vcf2", help="Second VCF file. Supports .vcf and .vcf.gz")

    parser.add_argument(
        "-o",
        "--out-prefix",
        default="vcf_compare",
        help="Output prefix. Default: vcf_compare",
    )

    args = parser.parse_args()

    variants1, stats1 = parse_vcf(Path(args.vcf1))
    variants2, stats2 = parse_vcf(Path(args.vcf2))

    keys1 = set(variants1)
    keys2 = set(variants2)

    shared_keys = keys1 & keys2
    unique1_keys = keys1 - keys2
    unique2_keys = keys2 - keys1

    unique1_variants = {key: variants1[key] for key in unique1_keys}
    unique2_variants = {key: variants2[key] for key in unique2_keys}

    out_prefix = args.out_prefix

    write_variant_table(f"{out_prefix}.unique_to_vcf1.tsv", unique1_variants)

    write_variant_table(f"{out_prefix}.unique_to_vcf2.tsv", unique2_variants)

    write_stats(
        f"{out_prefix}.stats.tsv",
        stats1,
        stats2,
        shared_count=len(shared_keys),
        unique1_count=len(unique1_keys),
        unique2_count=len(unique2_keys),
    )

    print("VCF comparison complete")
    print()
    print(f"VCF1: {args.vcf1}")
    print(f"VCF2: {args.vcf2}")
    print()
    print(f"VCF1 variants:        {len(keys1)}")
    print(f"VCF2 variants:        {len(keys2)}")
    print(f"Shared variants:      {len(shared_keys)}")
    print(f"Unique to VCF1:       {len(unique1_keys)}")
    print(f"Unique to VCF2:       {len(unique2_keys)}")
    print()
    print("Output files:")
    print(f"  {out_prefix}.stats.tsv")
    print(f"  {out_prefix}.unique_to_vcf1.tsv")
    print(f"  {out_prefix}.unique_to_vcf2.tsv")


if __name__ == "__main__":
    main()
