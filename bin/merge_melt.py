#!/usr/bin/env python3
"""Merge MELT mobile-element VCF files."""

import re
from argparse import ArgumentParser
from pathlib import Path

import vcf2 as vcf
from vcf_pipeline_utils import leading_float


def _natural_key(value):
    """Return a sort key close to GNU sort -V for chromosome strings."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def merge_melt(vcf_header, sample_id, out_file):
    """Merge MELT ALU, LINE1, and SVA VCF files into one Scout-style VCF."""
    qcval = {
        "lc": "Low Complex Region",
        "s25": "Greater than 25.0% of samples do not have data",
        "hDP": '"More than the expected number of discordant pairs at this site are also split',
        "PASS": "PASS",
        "rSD": "Ratio of LP to RP is greater than 2.0 standard deviations",
    }
    assess = {
        "0": "No overlapping reads at site",
        "1": "Imprecise breakpoint due to greater than expected distance between evidence",
        "2": "discordant pair evidence only -- No split read information",
        "3": "left side TSD evidence only",
        "4": "right side TSD evidence only",
        "5": "TSD decided with split reads -> highest possible quality",
    }
    mod_files = []
    for path in ("ALU.final_comp.vcf", "LINE1.final_comp.vcf", "SVA.final_comp.vcf"):
        if not Path(path).is_file() or Path(path).stat().st_size < 1:
            continue
        out_path = f"{path}mod"
        mod_files.append(out_path)
        with vcf.VCFReader(path) as reader, open(out_path, "w") as out_fh:
            inserted = False
            for line in reader.meta_header_lines:
                line = line.replace("'\\|'", "pipe-sign")
                if line.startswith("##INFO") and not inserted:
                    print(
                        '##INFO=<ID=END,Number=.,Type=Integer,Description="END position set to start position for insertions">',
                        file=out_fh,
                    )
                    print(
                        '##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description="Custom annotations for scout">',
                        file=out_fh,
                    )
                    inserted = True
                print(line, file=out_fh)
            print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}", file=out_fh)
            for var in reader:
                if "ac0" in var["FILTER"]:
                    continue
                alt_parts = var["ALT"].split(":")
                alt = (alt_parts[2] if len(alt_parts) > 2 else var["ALT"]).replace(">", "")
                event = str(var["INFO"].get("INTERNAL", "")).replace(",", "_")
                meinfo = str(var["INFO"].get("MEINFO", "")).replace(",", " ")
                custom = [
                    f"SCOUT_CUSTOM=Repeat Element|{alt}",
                    f"Subtype|{meinfo}",
                    f"Consequence|{event}",
                    f"Target Site Duplication|{var['INFO'].get('TSD', '')}",
                    f"Filter|{qcval.get(var['FILTER'], '')}",
                    f"Assess|{assess.get(str(var['INFO'].get('ASSESS', '')), '')}",
                ]
                gt = var["GT"][0]
                depth = leading_float(gt.get("DP"))
                alt_depth = leading_float(gt.get("AD"))
                vaf_value = f"{alt_depth / depth:.4f}" if depth else "0.0000"
                row = [
                    var["CHROM"],
                    var["POS"],
                    var["ID"],
                    var["REF"],
                    var["ALT"],
                    var["QUAL"],
                    "PASS",
                    f"SVTYPE=INS;END={var['POS']};{','.join(custom)}",
                    "GT:VAF:VD:DP",
                    f"{gt.get('GT', '')}:{vaf_value}:{gt.get('AD', '')}:{gt.get('DP', '')}",
                ]
                print("\t".join(row), file=out_fh)

    if not mod_files:
        with open(vcf_header) as header_fh, open(out_file, "w") as out_fh:
            for line in header_fh:
                out_fh.write(line.replace("FORMAT", f"FORMAT\t{sample_id}") if line.startswith("#CHROM") else line)
        return out_file

    lines = []
    header = []
    for idx, path in enumerate(mod_files):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    if idx == 0:
                        header.append(line)
                else:
                    lines.append(line)
    lines.sort(key=lambda row: (_natural_key(row.split("\t", 1)[0]), int(row.split("\t", 2)[1])))
    with open(out_file, "w") as out_fh:
        out_fh.writelines(header + lines)
    return out_file


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-header", required=True, help="Header VCF used when no records are present.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--out", required=True, help="Output merged VCF.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    merge_melt(args.vcf_header, args.sample_id, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
