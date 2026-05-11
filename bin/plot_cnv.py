#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

###############################################################################
# Custom CNV plotting:
#   - Top panel: denoised copy ratios from gatk or cnvkit
#   - Bottom panel: B-allele fractions (BAF)
#   - Optional segmentation overlays
#
# Input:
#   1. denoisedCR.tsv/cnvkit.cnr  (from GATK DenoiseReadCounts/cnvkit batch)
#   2. hets.tsv                   (heterozygous allelic counts)
#
# Example:
#
# python plot_cnv.py \
#     --cr sample.denoisedCR.tsv \
#     --hets sample.hets.tsv \
#     --out sample_cnv.png
#
###############################################################################


def normalize_chr(x):
    x = str(x)

    if x.startswith("chr"):
        x = x[3:]

    return x


def detect_column(df, candidates):
    """
    Find first matching column name from candidate list.
    """
    cols_lower = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]

    raise ValueError(
        f"Could not find any of columns: {candidates}\n"
        f"Available columns:\n{list(df.columns)}"
    )


def compute_baf(df, ref_col, alt_col):
    total = df[ref_col] + df[alt_col]
    baf = df[alt_col] / total
    return baf


def build_genome_positions(df, chrom_col, pos_col, chrom_order):
    """
    Create cumulative genomic coordinates.
    """

    chrom_sizes = df.groupby(chrom_col)[pos_col].max().reindex(chrom_order).fillna(0)

    offsets = {}
    running = 0

    for chrom in chrom_order:
        offsets[chrom] = running
        running += chrom_sizes[chrom]

    df["genome_pos"] = df[pos_col] + df[chrom_col].map(offsets)

    return df, offsets, chrom_sizes


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-type",
        choices=["gatk", "cnvkit"],
        default="gatk",
        help=(
            "Input format:\n"
            "  gatk   = GATK denoisedCR.tsv\n"
            "  cnvkit = CNVkit .cnr file"
        ),
    )

    parser.add_argument("--cr", required=True, help="denoisedCR.tsv")

    parser.add_argument("--hets", required=True, help="hets.tsv")

    parser.add_argument("--out", required=True, help="output image")

    parser.add_argument(
        "--sample-id", required=True, help="Sample ID to use as plot title"
    )

    parser.add_argument(
        "--segments", default=None, help="Optional CNVkit .cns segmentation file"
    )

    parser.add_argument("--point-size", type=float, default=3)

    parser.add_argument("--alpha", type=float, default=0.6)

    parser.add_argument(
        "--min-allele-depth",
        type=int,
        default=100,
        help="Minimum total allele depth for BAF plotting"
    )

    args = parser.parse_args()

    ###########################################################################
    # LOAD COPY RATIOS
    ###########################################################################

    cr = pd.read_csv(args.cr, sep="\t", comment="@")

    ###########################################################################
    # GATK FORMAT
    ###########################################################################

    if args.input_type == "gatk":

        chrom_col = detect_column(cr, ["CONTIG", "Chromosome", "chr"])

        start_col = detect_column(cr, ["START", "POSITION", "POS"])

        cr_col = detect_column(
            cr, ["LOG2_COPY_RATIO", "LOG2_COPY_RATIO_POSTERIOR_50", "COPY_RATIO"]
        )

    ###########################################################################
    # CNVKIT FORMAT (.cnr)
    ###########################################################################

    elif args.input_type == "cnvkit":

        chrom_col = detect_column(cr, ["chromosome"])

        start_col = detect_column(cr, ["start"])

        cr_col = detect_column(cr, ["log2"])

    ###########################################################################
    # CHROMOSOME ORDER
    ###########################################################################

    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]

    ###########################################################################
    # STANDARDIZE COPY-RATIO CHROMOSOME NAMES
    ###########################################################################

    cr[chrom_col] = cr[chrom_col].map(normalize_chr)

    cr = cr[cr[chrom_col].isin(chrom_order)].copy()

    ###########################################################################
    # LOAD OPTIONAL SEGMENTS
    ###########################################################################

    segments = None

    if args.segments is not None:

        segments = pd.read_csv(args.segments, sep="\t", comment="@")

        seg_chrom_col = detect_column(segments, ["chromosome", "CONTIG", "chr"])

        seg_start_col = detect_column(segments, ["start", "START"])

        seg_end_col = detect_column(segments, ["end", "END"])

        seg_log2_col = detect_column(segments, ["log2", "LOG2_COPY_RATIO"])

        # normalize chromosome naming
        segments[seg_chrom_col] = segments[seg_chrom_col].map(normalize_chr)

        # keep canonical chromosomes only
        segments = segments[segments[seg_chrom_col].isin(chrom_order)].copy()

    ###########################################################################
    # LOAD HETS
    ###########################################################################

    hets = pd.read_csv(
        args.hets,
        sep="\t",
        comment="@",
        dtype=str,
        low_memory=False
    )

    h_chrom_col = detect_column(hets, ["CONTIG", "Chromosome", "chr"])

    h_pos_col = detect_column(hets, ["POSITION", "POS", "START"])

    ref_col = detect_column(hets, ["REF_COUNT", "REFERENCE_COUNT"])

    alt_col = detect_column(hets, ["ALT_COUNT", "ALTERNATE_COUNT"])

    ###########################################################################
    # CONVERT NUMERIC COLUMNS
    ###########################################################################

    hets[h_pos_col] = pd.to_numeric(hets[h_pos_col])

    hets[ref_col] = pd.to_numeric(hets[ref_col])

    hets[alt_col] = pd.to_numeric(hets[alt_col])

    ###########################################################################
    # STANDARDIZE HETS CHROMOSOME NAMES
    ###########################################################################

    hets[h_chrom_col] = hets[h_chrom_col].map(normalize_chr)

    hets = hets[hets[h_chrom_col].isin(chrom_order)].copy()

    ###########################################################################
    # GENOME COORDINATES
    ###########################################################################

    cr, offsets, chrom_sizes = build_genome_positions(
        cr, chrom_col, start_col, chrom_order
    )

    hets["genome_pos"] = hets[h_pos_col] + hets[h_chrom_col].map(offsets)

    ###########################################################################
    # SEGMENT GENOMIC COORDINATES
    ###########################################################################

    if segments is not None:

        segments["genome_start"] = segments[seg_start_col] + segments[
            seg_chrom_col
        ].map(offsets)

        segments["genome_end"] = segments[seg_end_col] + segments[seg_chrom_col].map(
            offsets
        )

        # convert log2 -> linear copy ratio
        segments["copy_ratio_linear"] = 2 ** segments[seg_log2_col]

    ###########################################################################
    # COMPUTE BAF
    ###########################################################################

    ###########################################################################
    # COMPUTE TOTAL DEPTH
    ###########################################################################

    hets["TOTAL_DEPTH"] = (
        hets[ref_col] +
        hets[alt_col]
    )

    ###########################################################################
    # FILTER LOW-DEPTH LOCI
    ###########################################################################

    hets = hets[
        hets["TOTAL_DEPTH"] >= args.min_allele_depth
    ].copy()

    ###########################################################################
    # COMPUTE BAF
    ###########################################################################

    hets["BAF"] = compute_baf(
        hets,
        ref_col,
        alt_col
    )

    ###########################################################################
    # PLOT
    ###########################################################################

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(20, 10), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    ###########################################################################
    # CONVERT LOG2 -> LINEAR COPY RATIO
    ###########################################################################

    cr["COPY_RATIO_LINEAR"] = 2 ** (cr[cr_col])

    ###########################################################################
    # CAP EXTREME AMPLIFICATIONS
    ###########################################################################

    MAX_COPY_RATIO = 4

    cr["COPY_RATIO_PLOT"] = (
        cr["COPY_RATIO_LINEAR"]
        .clip(upper=MAX_COPY_RATIO)
    )

    ###########################################################################
    # IDENTIFY CAPPED POINTS
    ###########################################################################

    high_amp = (
        cr["COPY_RATIO_LINEAR"] > MAX_COPY_RATIO
    )

    normal_amp = ~high_amp

    ###########################################################################
    # COPY RATIO PANEL
    ###########################################################################

    # normal points
    ax1.scatter(
        cr.loc[normal_amp, "genome_pos"],
        cr.loc[normal_amp, "COPY_RATIO_PLOT"],
        s=args.point_size,
        alpha=args.alpha,
        rasterized=True
    )

    # capped amplification points
    ax1.scatter(
        cr.loc[high_amp, "genome_pos"],
        cr.loc[high_amp, "COPY_RATIO_PLOT"],
        s=args.point_size * 4,
        marker="^",
        alpha=args.alpha,
        rasterized=True
    )

    # Neutral baseline
    ax1.axhline(1.0, linestyle="--")

    ax1.set_ylabel("Copy Ratio")
    ax1.set_title(args.sample_id)

    # Fixed y-axis range
    ax1.set_ylim(0, 4)

    ###########################################################################
    # OPTIONAL SEGMENT OVERLAY
    ###########################################################################

    if segments is not None:

        for _, row in segments.iterrows():

            ax1.plot(
                [row["genome_start"], row["genome_end"]],
                [row["copy_ratio_linear"], row["copy_ratio_linear"]],
                linewidth=2,
                color="black",
            )

    ###########################################################################
    # BAF PANEL
    ###########################################################################

    ax2.scatter(
        hets["genome_pos"],
        hets["BAF"],
        s=args.point_size,
        alpha=args.alpha,
        rasterized=True,
    )

    ax2.axhline(0.5, linestyle="--")

    ax2.set_ylabel("BAF")
    ax2.set_ylim(0, 1)

    ###########################################################################
    # CHROMOSOME LABELS
    ###########################################################################

    tick_positions = []
    tick_labels = []

    running = 0

    for chrom in chrom_order:

        size = chrom_sizes.get(chrom, 0)

        if size == 0:
            continue

        midpoint = running + size / 2

        tick_positions.append(midpoint)
        tick_labels.append(chrom)

        ax1.axvline(running, linewidth=0.5)
        ax2.axvline(running, linewidth=0.5)

        running += size

    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(tick_labels)

    ax2.set_xlabel("Chromosome")

    ###########################################################################
    # CLEANUP
    ###########################################################################

    plt.tight_layout()

    plt.savefig(args.out, dpi=300, bbox_inches="tight")

    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
