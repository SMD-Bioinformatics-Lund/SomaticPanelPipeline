"""CNV and coverage utility commands."""

from __future__ import annotations

import json
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import vcf2 as vcf
from vcf_pipeline_utils import leading_float, open_output


def cnvkit_to_hrd(input_file, sample_id, ploidy, out_file):
    """Convert CNVkit call table to HRD input format."""
    with open(input_file) as in_fh, open_output(out_file) as out_fh:
        suffix = "\tploidy" if ploidy else ""
        out_fh.write(f"SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn{suffix}\n")
        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith("chromosome") or not line:
                continue
            fields = line.split("\t")
            if fields[0] == "Y" and leading_float(fields[8]) == 0:
                continue
            a_cn = str(leading_float(fields[8]) - 1) if fields[10] == "" or fields[9] == "" else fields[10]
            b_cn = "1" if fields[10] == "" or fields[9] == "" else fields[9]
            row = [sample_id, f"chr{fields[0]}", fields[1], fields[2], fields[8], a_cn, b_cn]
            if ploidy:
                row.append(ploidy)
            out_fh.write("\t".join(row) + "\n")


def filter_cnvkit(tsv_file, coverage, distance, out_file):
    """Filter CNVkit TSV calls."""
    cov = leading_float(coverage)
    max_dist = leading_float(distance)
    with open(tsv_file) as in_fh, open_output(out_file) as out_fh:
        for idx, line in enumerate(in_fh, start=1):
            fields = line.split("\t")
            if idx == 1:
                out_fh.write(line)
                continue
            if fields[3] == "-":
                continue
            if leading_float(fields[6]) == 2:
                continue
            if leading_float(fields[2]) - leading_float(fields[1]) > max_dist:
                continue
            if leading_float(fields[6]) == 1 and cov and leading_float(fields[9]) / cov <= 0.10:
                continue
            out_fh.write("\t".join(fields))


def gatk_segments(cr_file, purity, out_file):
    """Convert GATK segment table to fold-change table."""
    with open(cr_file) as in_fh, open_output(out_file) as out_fh:
        for line in in_fh:
            if line.startswith("@") or line.startswith("C"):
                continue
            fields = line.rstrip("\n").split("\t")
            foldchange = (2 ** leading_float(fields[4])) / leading_float(purity)
            out_fh.write(f"{fields[0]}\t{fields[4]}\t{foldchange}\t{fields[5]}\n")


def generate_gens_data(cnr_file, vcf_file, sample_id, out_prefix):
    """Generate compressed CNVkit coverage and BAF BED tracks."""
    log2_data = []
    with open(cnr_file) as cnr:
        next(cnr, None)
        for line in cnr:
            chrom, start, end, gene, _depth, log2, *_rest = line.rstrip("\n").split("\t")
            if gene == "Antitarget":
                continue
            midpoint = int(start) + int((int(end) - int(start)) / 2)
            log2_data.append(f"{chrom}\t{midpoint - 1}\t{midpoint}\t{log2}")

    baf_data = []
    with vcf.VCFReader(vcf_file) as reader:
        for var in reader:
            for gt in var["GT"]:
                if gt["_sample_id"] != sample_id:
                    continue
                depth = leading_float(gt.get("DP"))
                if depth < 100 or not gt.get("AO"):
                    break
                vaf_value = leading_float(gt.get("AO")) / depth
                baf_data.append(f"{var['CHROM']}\t{int(var['POS']) - 1}\t{var['POS']}\t{vaf_value}")

    prefix = out_prefix or sample_id
    cov_bed = f"{prefix}.cov.bed"
    baf_bed = f"{prefix}.baf.bed"
    for path, data in ((cov_bed, log2_data), (baf_bed, baf_data)):
        with open(path, "w") as out_fh:
            for resolution in "oabcd":
                for row in data:
                    out_fh.write(f"{resolution}_{row}\n")
        subprocess.run(["bgzip", "-f", path], check=True)
        subprocess.run(["tabix", f"{path}.gz"], check=True)


def panel_depth(bam_file, bed_file, out_file, cutoff=500):
    """Report contiguous low-depth panel regions from sambamba depth output."""
    if not Path(bam_file).is_file():
        raise FileNotFoundError(f"No file {bam_file}")
    if not Path(bed_file).is_file():
        raise FileNotFoundError(f"No file {bed_file}")
    cmd = ["sambamba", "depth", "base", bam_file, "-L", bed_file]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc, open_output(out_file) as out_fh:
        start_pos = start_chr = last_low_pos = last_low_chr = None
        low_cov_sum = 0.0
        for line in proc.stdout:
            chrom, pos_text, depth_text, *_rest = line.split("\t")
            pos = int(pos_text)
            depth = leading_float(depth_text)
            if depth < cutoff:
                if last_low_chr == chrom and last_low_pos == pos - 1:
                    low_cov_sum += depth
                else:
                    if start_pos is not None and start_chr is not None:
                        out_fh.write(
                            f"{start_chr}\t{start_pos}\t{last_low_pos}\t{low_cov_sum / (last_low_pos - start_pos + 1)}\n"
                        )
                    start_chr = chrom
                    start_pos = pos
                    low_cov_sum = depth
                last_low_chr = chrom
                last_low_pos = pos
            elif start_chr is not None and start_pos is not None:
                out_fh.write(
                    f"{start_chr}\t{start_pos}\t{last_low_pos}\t{low_cov_sum / (last_low_pos - start_pos + 1)}\n"
                )
                start_chr = start_pos = None
        if proc.wait() != 0:
            raise RuntimeError("sambamba depth base failed")


def qc_sentieon(sample_id, assay_type, out_file):
    """Parse Sentieon QC metric files in the current directory."""
    results = {}
    pct_above = {}
    if assay_type == "panel":
        pct_above, median = panel_coverage_calc()
        results["median_cov"] = median
        hs_file = "hs_metrics.txt"
        for row in metric_value_rows(hs_file):
            results["pct_on_target"] = row[18]
            results["fold_enrichment"] = row[25]
            results["mean_coverage"] = row[22]
            results["fold_80"] = row[32]
            results["at_drop"] = row[48]
            results["gc_drop"] = row[49]
    elif assay_type == "wgs":
        pct_above = parse_wgs_metrics(results)
        parse_gc_summary(results)
    else:
        raise ValueError("assay type must be panel or wgs")

    for row in metric_value_rows("is_metrics.txt"):
        results["ins_size"] = f"{leading_float(row[4]):.0f}"
        results["ins_size_dev"] = f"{leading_float(row[5]):.0f}"
    for row in metric_value_rows("dedup_metrics.txt"):
        results["dup_reads"] = row[6]
        results["num_reads"] = row[2]
        results["dup_pct"] = row[8]
        results["mapped_reads"] = leading_float(row[2]) - leading_float(row[4])
    for row in metric_value_rows("aln_metrics.txt"):
        results["pf_mismatch_rate"] = row[12]
        results["pf_error_rate"] = row[13]

    results["pct_above_x"] = pct_above
    results["sample_id"] = sample_id
    with open_output(out_file) as out_fh:
        json.dump(results, out_fh, indent=3)
        out_fh.write("\n")


def metric_value_rows(path):
    """Yield Sentieon metric value rows following command/header lines."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#SentieonCommandLine"):
                next(fh, None)
                yield next(fh).rstrip("\n").split("\t")


def panel_coverage_calc():
    """Calculate panel coverage summary values."""
    cov = []
    with open("cov_metrics.txt") as fh:
        for line in fh:
            if not line.startswith("Locus"):
                vals = next(fh, "").split("\t")
                if len(vals) > 1:
                    cov.append(leading_float(vals[1]))
    sorted_cov = sorted(cov)
    mid = len(sorted_cov) // 2
    median = sorted_cov[mid] if len(sorted_cov) % 2 else sorted_cov[mid - 1]
    pct = {}
    with open("cov_metrics.txt.sample_summary") as fh:
        for line in fh:
            if line.startswith("sample_id"):
                vals = next(fh).rstrip("\n").split("\t")
                for key, idx in (("1", 6), ("10", 7), ("30", 8), ("100", 9), ("250", 10), ("500", 11)):
                    pct[key] = vals[idx]
    return pct, median


def parse_wgs_metrics(results):
    """Parse WGS metrics and return pct-above thresholds."""
    pct = {}
    quartiles = {}
    total = 0
    pct_obs = {}
    with open("wgs_metrics.txt") as fh:
        for line in fh:
            if line.startswith("#SentieonCommandLine"):
                next(fh, None)
                row = next(fh).split("\t")
                results["median_cov"] = row[3]
                results["sd_coverage"] = row[2]
                results["mean_coverage"] = row[1]
                pct_obs = {
                    "R_25": (leading_float(row[0]) + 1) / 4,
                    "R_50": (leading_float(row[0]) + 1) / 2,
                    "R_75": ((leading_float(row[0]) + 1) * 3) / 4,
                }
                for key, idx in zip(
                    ("1", "5", "10", "15", "20", "25", "30", "40", "50", "60", "70", "80", "90", "100"),
                    range(12, 26),
                    strict=False,
                ):
                    pct[key] = row[idx]
            elif line[:1].isdigit():
                depth, count = line.split()[:2]
                total += leading_float(count)
                for key, obs in pct_obs.items():
                    if total >= obs and key not in quartiles:
                        quartiles[key] = leading_float(depth)
    results["iqr"] = quartiles.get("R_75", 0) - quartiles.get("R_25", 0)
    return pct


def parse_gc_summary(results):
    """Parse GC summary file into results."""
    for row in metric_value_rows("gc_summary.txt"):
        results["at_drop"] = row[5]
        results["gc_drop"] = row[6]


def build_parser():
    parser = ArgumentParser(description="Report contiguous low-depth panel regions.")
    parser.add_argument("--bam", required=True, help="Input BAM.")
    parser.add_argument("--bed", required=True, help="Panel BED.")
    parser.add_argument("--cutoff", type=float, default=500, help="Low-depth cutoff.")
    parser.add_argument("--out", required=True, help="Output BED.")
    return parser

def main(argv=None):
    args = build_parser().parse_args(argv)
    panel_depth(args.bam, args.bed, args.out, args.cutoff)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
