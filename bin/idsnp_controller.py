"""ID-SNP concordance utilities."""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

import vcf2 as vcf
from vcf_pipeline_utils import leading_float


def _dp4_from_info(info_text):
    """Return ref and alt depths from a VCF INFO DP4 field."""
    for item in info_text.split(";"):
        if item.startswith("DP4="):
            vals = [leading_float(value) for value in item[4:].split(",")]
            while len(vals) < 4:
                vals.append(0)
            return vals[0] + vals[1], vals[2] + vals[3]
    return 0, 0


def _read_idsnp_vcf(path, sample_label):
    """Read one ID-SNP VCF into a marker dictionary."""
    data = {}
    with vcf.open_vcf(path) as fh:
        labels = {}
        for line in fh:
            line = line.rstrip("\n")
            cols = line.split("\t")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                labels = {value: idx for idx, value in enumerate(cols)}
                continue
            identifier = f"{cols[labels['#CHROM']]}_{cols[labels['POS']]}_{cols[labels['ID']]}"
            gt = cols[labels["FORMAT"] + 1].split(":")[0]
            dp_ref, dp_alt = _dp4_from_info(cols[labels["INFO"]])
            data[identifier] = {
                sample_label: {
                    "genotype": gt,
                    "quality": cols[labels["QUAL"]],
                    "dp": dp_ref + dp_alt,
                    "dp_ref": dp_ref,
                    "dp_alt": dp_alt,
                }
            }
    return data


def _merge_marker_data(target, source):
    """Merge marker dictionaries in place."""
    for key, value in source.items():
        target.setdefault(key, {}).update(value)


def compare_idsnps(vcf_sample, vcf_control, sample, control, rs_bed, out_prefix=None):
    """Compare sample and control ID-SNP genotypes."""
    if not Path(vcf_sample).exists():
        raise FileNotFoundError("SAMPLE VCF file not found")
    if not Path(vcf_control).exists():
        raise FileNotFoundError(f"Control VCF file '{vcf_control}' does not exist")

    marker_count = sum(1 for _line in open(rs_bed))
    data = {}
    _merge_marker_data(data, _read_idsnp_vcf(vcf_sample, sample))
    _merge_marker_data(data, _read_idsnp_vcf(vcf_control, control))

    prefix = out_prefix or f"s{sample}_c{control}"
    csv_path = f"{prefix}.csv"
    json_path = f"{prefix}.json"

    total = match = mismatch = 0
    with open(csv_path, "w") as out_fh:
        out_fh.write(
            "VARIANT\tgt-SAMPLE\tQ-SAMPLE\tDP-SAMPLE[ref/alt]\tgt-CONTROL\tQ-CONTROL\tDP-CONTROL[ref/alt]\tCONCORDANCE\n"
        )
        for marker in sorted(data):
            total += 1
            out_fh.write(f"{marker}\t")
            sample_data = data[marker].get(sample)
            control_data = data[marker].get(control)
            if sample_data:
                out_fh.write(
                    f"{sample_data['genotype']}\t{sample_data['quality']}\t"
                    f"{sample_data['dp']}[{sample_data['dp_ref']}/{sample_data['dp_alt']}]"
                )
            else:
                out_fh.write("NA\tNA\tNA")
            if control_data:
                out_fh.write(
                    f"\t{control_data['genotype']}\t{control_data['quality']}\t"
                    f"{control_data['dp']}[{control_data['dp_ref']}/{control_data['dp_alt']}]"
                )
            else:
                out_fh.write("\tNA\tNA\tNA")
            if sample_data and control_data and sample_data["genotype"] == control_data["genotype"]:
                out_fh.write("\tMATCH\n")
                match += 1
            else:
                out_fh.write("\tMISMATCH\n")
                mismatch += 1

        missing = marker_count - total
        out_fh.write(
            f"\n# missing data: {missing} calls out of {marker_count} ({(missing / marker_count) * 100:.2f})\n"
        )
        out_fh.write(f"# mismatches: {mismatch}\n")
        out_fh.write(f"# found matches: {match} out of available markers({total}) ({(match / total) * 100:.2f})\n")
        out_fh.write(
            f"# total matches: {match} out of all markers({marker_count}) ({(match / marker_count) * 100:.2f})\n"
        )

    result = {
        "is_paired_sample": True,
        "pct_matching_snps": float(f"{(match / marker_count) * 100:.2f}"),
        "nbr_non_matching_snps": int(mismatch),
        "nbr_matching_snps": int(match),
        "nbr_missing_snps": int(missing),
        "total_nbr_id_snps": int(total),
    }
    with open(json_path, "w") as out_fh:
        json.dump(result, out_fh, indent=3)
        out_fh.write("\n")
    return csv_path, json_path


def build_parser():
    """Build the ID-SNP parser."""
    parser = ArgumentParser(description="Compare ID-SNP genotypes between sample and control VCFs.")
    parser.add_argument("--vcf-sample", required=True, help="Sample VCF.")
    parser.add_argument("--vcf-control", required=True, help="Control VCF.")
    parser.add_argument("--sample", required=True, help="Sample ID.")
    parser.add_argument("--control", required=True, help="Control ID.")
    parser.add_argument("--rs-bed", required=True, help="BED file of ID-SNP markers.")
    parser.add_argument("--out-prefix", help="Output prefix. Defaults to s<SAMPLE>_c<CONTROL>.")
    return parser


def main(argv=None):
    """Run the ID-SNP command."""
    args = build_parser().parse_args(argv)
    _csv, json_path = compare_idsnps(
        args.vcf_sample, args.vcf_control, args.sample, args.control, args.rs_bed, args.out_prefix
    )
    print(f"Data has been written to {json_path}")
    return 0
