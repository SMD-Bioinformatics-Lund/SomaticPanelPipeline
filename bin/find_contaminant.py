#!/usr/bin/env python3
"""Estimate sample contamination from common gnomAD SNV VAF distributions."""

import os
import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf2  # noqa: E402
from vcf_pipeline_utils import leading_float, open_output  # noqa: E402


EMPTY_PNG = (
    b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
    b"\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc```\x00\x00"
    b"\x00\x04\x00\x01\xf6\x178U\x00\x00\x00\x00IEND\xaeB`\x82"
)


def max_ampersand_float(value):
    """Return the maximum numeric value from an ampersand-delimited string."""
    maximum = 0.0
    for item in str(value or "").split("&"):
        number = leading_float(item)
        if number > maximum:
            maximum = number
    return maximum


def paired_sample(case_id, var):
    """Return the paired sample ID when the VCF has more than one sample."""
    normal_id = 0
    if len(var["GT"]) > 1:
        for gt in var["GT"]:
            if gt["_sample_id"] != case_id:
                normal_id = gt["_sample_id"]
    return normal_id


def mean(values):
    """Return the arithmetic mean of values."""
    return sum(values) / len(values)


def write_placeholder_png(path):
    """Write a minimal PNG file.

    Args:
        path: PNG path.
    """
    with open(path, "wb") as fh:
        fh.write(EMPTY_PNG)


class ContaminationEstimator:
    """Stateful contamination estimator matching the legacy window logic."""

    def __init__(self, case_id, low=0.0, high=0.3, ad_field="VD", binsize_cutoff=80):
        """Initialize estimator state.

        Args:
            case_id: Sample ID to score.
            low: Lower VAF detection limit.
            high: Upper VAF limit.
            ad_field: FORMAT field used for alternate depth.
            binsize_cutoff: Minimum peak bin size.
        """
        self.case_id = case_id
        self.low = low
        self.high = high
        self.ad_field = ad_field
        self.binsize_cutoff = binsize_cutoff
        self.dist_dec = {}

    def add_vaf(self, vaf, other_vaf, var_key):
        """Add a VAF to the distribution if it passes legacy checks.

        Args:
            vaf: VAF being evaluated.
            other_vaf: Paired sample VAF, or zero for unpaired mode.
            var_key: Tab-delimited variant key used in genotype output.
        """
        if not other_vaf:
            if self.low <= vaf <= self.high:
                bucket = self.dist_dec.setdefault(vaf, {"COUNT": 0, "VAR": []})
                bucket["COUNT"] += 1
                bucket["VAR"].append(var_key)
        elif vaf <= self.low and self.low <= other_vaf <= self.high:
            bucket = self.dist_dec.setdefault(other_vaf, {"COUNT": 0, "VAR": []})
            bucket["COUNT"] += 1
            bucket["VAR"].append(var_key)

    def distribution(self):
        """Build sliding-window VAF distribution bins.

        Returns:
            Tuple of distribution dict, number of bins, and total variants in
            all non-empty bins.
        """
        distri = {}
        inc = 0.005
        window = 0.005
        num_vars_bin = 0
        num_bins = 0
        current = self.low
        while current <= self.high + 1e-12:
            total = 0
            afs = []
            variants = []
            for af in sorted(self.dist_dec):
                if current <= af <= current + window:
                    total += self.dist_dec[af]["COUNT"]
                    afs.append(af)
                    variants.extend(self.dist_dec[af]["VAR"])
            if total:
                num_bins += 1
                distri[num_bins] = {
                    "AFspan": f"{current}-{current + window}",
                    "COUNT": total,
                    "MEAN": mean(afs),
                    "VARS": variants,
                }
                num_vars_bin += total
            current += inc
        return distri, num_bins, num_vars_bin

    def heterozygous_peak(self, distri, num_bins, num_vars_bin, output_id):
        """Find the likely heterozygous contaminant peak.

        Args:
            distri: Distribution bins.
            num_bins: Number of non-empty bins.
            num_vars_bin: Total variants across non-empty bins.
            output_id: Prefix for side-effect files.

        Returns:
            Tuple of peak AF and peak bin index.
        """
        dist_file = f"{output_id}.dist.txt"
        mean_bincount = num_vars_bin / num_bins
        highpoint = 0
        af_at_highpoint = 0
        bin_at_highpoint = 0
        with open(dist_file, "w") as out_fh:
            for af_count in sorted(distri):
                count = distri[af_count]["COUNT"]
                if count / mean_bincount >= 3.5 and count > highpoint and count > self.binsize_cutoff:
                    highpoint = count
                    af_at_highpoint = distri[af_count]["MEAN"]
                    bin_at_highpoint = af_count
                out_fh.write(f"{af_count}\t{count}\t{distri[af_count]['MEAN']}\n")
        write_placeholder_png(f"{output_id}.png")
        return af_at_highpoint, bin_at_highpoint

    @staticmethod
    def homozygous_peak(distri, first_peak, af_at_hetero, num_bins):
        """Find a possible homozygous peak around twice the heterozygous AF."""
        homo_highpoint = 0
        af_at_homo = 0
        bin_at_homo_highpoint = 0
        for idx in range(int(first_peak) + 1, num_bins):
            current = distri.get(idx, {}).get("COUNT", 0)
            next_count = distri.get(idx + 1, {}).get("COUNT", 0)
            if current < next_count and current > homo_highpoint:
                homo_highpoint = next_count
                af_at_homo = distri.get(idx + 1, {}).get("MEAN", 0)
                bin_at_homo_highpoint = idx + 1
        if af_at_hetero and 1.8 <= af_at_homo / af_at_hetero <= 2.2:
            return bin_at_homo_highpoint, homo_highpoint, af_at_homo
        return "cannot find homo", None, None

    @staticmethod
    def write_genotypes(output_id, variants, genotype_type):
        """Append genotype suggestions for peak variants."""
        with open(f"{output_id}.genotypes.txt", "a") as out_fh:
            for var in variants:
                chrom, pos, ref, alt = var.split("\t")
                if genotype_type == "hetero":
                    out_fh.write(f"{chrom}\t{pos}\t{ref}/{alt}\n")
                elif genotype_type == "homo":
                    out_fh.write(f"{chrom}\t{pos}\t{alt}/{alt}\n")


def genotype_vaf(gt, ad_field, current_tumor_vd=None):
    """Calculate VAF and depth from a genotype record.

    Args:
        gt: Parsed genotype dictionary.
        ad_field: FORMAT field used for allele depth.
        current_tumor_vd: Tumor AD value used to match a legacy normal branch.

    Returns:
        Tuple of VAF, allele depth field, and total depth.
    """
    allele_depth = gt.get(ad_field, "")
    depth = leading_float(gt.get("DP"))
    if ad_field == "VD":
        return leading_float(gt.get("VAF")), allele_depth, depth
    depth_values = str(current_tumor_vd if current_tumor_vd is not None else allele_depth).split(",")
    minimum = min(leading_float(item) for item in depth_values) if depth_values else 0.0
    return (0.0 if depth == 0 else minimum / depth), allele_depth, depth


def estimate_contamination(vcf_file, case_id, check_normal, detect_level, ad_field, high, binsize_cutoff, out_file):
    """Estimate contamination and write result and side files.

    Args:
        vcf_file: Input VCF path.
        case_id: Sample ID to evaluate.
        check_normal: Whether to output the paired normal score.
        detect_level: Lower detection level.
        ad_field: FORMAT field used for allele depth.
        high: Upper VAF limit.
        binsize_cutoff: Minimum bin count for peak calling.
        out_file: Output value path.
    """
    estimator = ContaminationEstimator(
        case_id=case_id,
        low=detect_level,
        high=high,
        ad_field=ad_field,
        binsize_cutoff=binsize_cutoff,
    )
    paired_id = 0
    messages = []

    with vcf2.VCFReader(file=vcf_file) as reader:
        for var in reader:
            if "X" in var["CHROM"] or "Y" in var["CHROM"]:
                continue
            csq = var["INFO"].get("CSQ") or []
            if not csq or not str(csq[0].get("gnomADg_AF", "")).strip():
                continue
            gnomad = max_ampersand_float(csq[0].get("gnomADg_AF"))
            if gnomad < 0.05:
                continue
            if len(var["REF"]) > 1 or len(var["ALT"]) > 1:
                continue

            var_key = f"{var['CHROM']}\t{var['POS']}\t{var['REF']}\t{var['ALT']}"
            paired_id = paired_sample(case_id, var)
            if check_normal and paired_id == 0:
                messages.append("no normal sample in vcf\n")

            tumor_vaf = tumor_depth = normal_vaf = 0.0
            tumor_vd = None
            for gt in var["GT"]:
                if gt["_sample_id"] == case_id:
                    tumor_vaf, tumor_vd, tumor_depth = genotype_vaf(gt, ad_field)
                else:
                    normal_vaf, _normal_vd, _normal_depth = genotype_vaf(gt, ad_field, current_tumor_vd=tumor_vd)

            if tumor_depth <= 5:
                continue
            if paired_id:
                if check_normal:
                    estimator.add_vaf(normal_vaf, 0, var_key)
                else:
                    estimator.add_vaf(normal_vaf, tumor_vaf, var_key)
            else:
                estimator.add_vaf(tumor_vaf, 0, var_key)

    output_id = paired_id if check_normal and paired_id else case_id
    distri, num_bins, num_vars_bin = estimator.distribution()
    with open_output(out_file) as out_fh:
        for message in messages:
            out_fh.write(message)
        if num_bins == 0:
            out_fh.write("0.0\n")
            Path(f"{output_id}.png").touch()
            return

        af_at_highpoint, bin_at_highpoint = estimator.heterozygous_peak(
            distri, num_bins, num_vars_bin, output_id
        )
        if af_at_highpoint:
            out_fh.write(f"{af_at_highpoint}\n")
            bin_at_homo, _homo_highpoint, _af_at_homo = estimator.homozygous_peak(
                distri, bin_at_highpoint, af_at_highpoint, num_bins
            )
            estimator.write_genotypes(output_id, distri[bin_at_highpoint]["VARS"], "hetero")
            if bin_at_homo != "cannot find homo":
                estimator.write_genotypes(output_id, distri[bin_at_homo]["VARS"], "homo")
        else:
            out_fh.write("0.0\n")


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Estimate contamination from VCF genotypes.")
    parser.add_argument("--vcf", required=True, help="Input annotated VCF file.")
    parser.add_argument("--case-id", required=True, help="Sample ID to score.")
    parser.add_argument("--normal", action="store_true", help="Score the paired normal sample.")
    parser.add_argument("--detect-level", type=float, default=0.0, help="Lower VAF detection level.")
    parser.add_argument("--ADfield-name", default="VD", help="FORMAT field used for allele depth.")
    parser.add_argument("--high", type=float, default=0.3, help="Upper VAF limit.")
    parser.add_argument("--binsize-cutoff", type=float, default=80, help="Minimum peak bin count.")
    parser.add_argument("--assay", help="Accepted for pipeline compatibility; ignored.")
    parser.add_argument("--out", required=True, help="Output contamination value file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    estimate_contamination(
        args.vcf,
        args.case_id,
        args.normal,
        args.detect_level,
        args.ADfield_name,
        args.high,
        args.binsize_cutoff,
        args.out,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
