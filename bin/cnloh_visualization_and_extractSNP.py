#!/usr/bin/env python3

"""
Description:
This script extracts informative SNPs from a posterior file based on LOH (Loss of Heterozygosity) regions provided in an input file.
These extracted SNPs are the ones used by PureCN for making LOH calls. For those extracted SNPs, a plot is generated for each chromosomes,
for the corresponding observed allele ratio (MAF)
The script also generates a bed files for GENS import for these CN_LOH calls.

Inputs:
1. --input_loh: A TSV file containing LOH calls with columns for chromosome, start, end, and LOH type.
2. --posterior_file: A TSV file containing SNP information with various attributes.
3. --output: A CSV file where the extracted informative SNPs will be written.

Outputs:
- A CSV file containing SNPs that meet the filtering criteria and fall within the specified CN_LOH regions.
- A TSV file for GENS visualization of the Copy Neutral LOH events (CN-LOH).
- A plot for each chromosome showing the B-allele frequency (BAF) for the extracted SNPs if CN_LOH calls are present.
- A JSON file for Coyote visualization of the Copy Neutral LOH events (CN-LOH).

Usage Example:
cnloh_visualization_and_extractSNP.py \
    --input_loh sample_loh_filtered.tsv \
    --posterior_file sample_snv_posteriors.tsv \
    --out_snp informative_snps.tsv \
    --out_gens sampleID_cnloh_gens.tsv
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extraction informative SNPs for final LOH calls"
    )

    parser.add_argument(
        '-i',
        '--input_loh',
        help='<sample>_loh_filtered.tsv is used to collect the chromosome and positions of the LOH calls',
        required=True
    )

    parser.add_argument(
        '-i2',
        '--posterior_file',
        help='<sample>_snv_posteriors.tsv is used to collect SNPs',
        required=True
    )

    parser.add_argument(
        '-o',
        '--out_snp',
        help='tsv file _snp.tsv containing the extracted informative SNPs used by PureCN for making LOH calls',
        required=False
    )

    parser.add_argument(
        '-o2',
        '--out_gens',
        help='bed file (9-column format) _cnloh_gens.bed for visualization of the Copy Neutral LOH events (CN-LOH)',
        required=True
    )

    return parser.parse_args()

args = parse_args()

# Generate the gens_tsv file with the CN-LOH segment for GENS visualization (violet segment)
cnloh_regions = {} # dictionary to collect LOH calls
annotations={}
with open(args.input_loh, "r") as loh_call, open(args.out_gens, "w") as gens_tsv:
    gens_tsv.write(f'chromosome\tstart\tend\tname\tcomments\tcolor\n')
    for line in loh_call:
        if not line.startswith("Sampleid"):
            field = line.split(",")
            chrom = field[1].strip()
            start = float(field[2].strip())
            end = float(field[3].strip())
            loh_type = field[7]
            snp = int(field[10].strip())
            delta_maf = float(field[14].strip())
            snp_per_mb = float(field[15].strip())
            purity=float(field[16].strip())
            purity_warning=field[17].strip()
            annotations.setdefault(chrom, []).append(
                (snp, delta_maf, snp_per_mb, purity, purity_warning)
            )
            if "COPY-NEUTRAL LOH" in loh_type:
                cnloh_regions.setdefault(chrom, []).append((start, end))

                comments = (
                    f"SNP: {snp};"
                    f"BAFexpected - BAF observed: {delta_maf};"
                    f"SNPs per Mb: {snp_per_mb};"
                    f"Purity: {purity};"
                    f"Purity accuracy warning: {purity_warning}"
                )
                gens_tsv.write(
                    "{chrom}\t{start}\t{end}\t"
                    "{loh_type}\t{comments}\trgb(127,0,255)\n".format(
                        chrom=chrom,
                        start=start,
                        end=end,
                        loh_type="CN-LOH",
                        comments=comments
                    )
                )

# Generate the file containing all informative SNPs for the CN-LOH calls used by PureCN
with open(args.posterior_file, "r") as posterior, open(args.out_snp, "w") as snps:
    for line in posterior:
        # Keep header
        if line.startswith("chr"):
            snps.write(line)
            continue

        field = line.split("\t")
        
        chrom = field[0].strip()
        position = int(field[1].strip())
        if chrom not in cnloh_regions:
            continue

        ml_somatic = field[25].strip()
        germline_contlow = float(field[23].strip())
        germline_conthigh = float(field[22].strip())
        germline_homozygous = float(field[24].strip())

        # --------------------------------------------------------------------------
        # Check SNP criteria
        # Likely germline, Likely heterozygous, Likely not cross-sample contamination
        # ---------------------------------------------------------------------------

        # remove somatic SNP (PureCN model's maximum-likelihood classification)
        if ml_somatic.lower() != "false":
            continue
        
        # probability that this variant is better explained by cross-sample contamination rather than a true germline genotype
            # CONTLOW for a contaminant peak at low VAF (we use 0.5 threshold as in PureCN own code in callLOH.R)
        if germline_contlow >= 0.5:
            continue
            # CONTHIGH for a contaminant peak at high VAF
        if germline_conthigh >= 0.5:
            continue
        
        # posterior probability of the SNP being germline homozygous
        if germline_homozygous >= 0.5:
            continue

        # Check whether SNP position is inside an LOH interval
        position_in_loh = any(
            start <= position <= end
            for start, end in cnloh_regions[chrom]
            )

        if not position_in_loh:
            continue

        # SNP passes all filters
        snps.write(line)


# MAKING THE PLOTs for BAF for the final CN LOH calls (one plot per chromosome)
# Read the snp table
df = pd.read_csv(args.out_snp, sep="\t")

# Get sample ID from output filename
sample_id = os.path.basename(args.out_snp).replace("_snp.tsv", "")

# Make one plot for each chromosome
for chromosome, data in df.groupby("chr"):

    plt.figure(figsize=(14, 6))

    #plt.scatter(data["start"], 
    #            data["ML.AR"], 
    #            marker="s", 
    #            facecolors="none",
    #            edgecolors="black",
    #            label="Exp allelic ratio (AR)")
    # plot BAF
    plt.scatter(data["start"], data["AR"], facecolors="hotpink", label="BAF")
    # plot BAF adjusted for mapping biais
    plt.scatter(data["start"], data["AR.ADJUSTED"], facecolors="darkgreen", label="BAF(mapping-biais adjusted)")

    # add a horizontal line at 0.5
    plt.axhline(
        0.5, 
        linestyle="--")

    # Labels and title
    plt.xlabel(
        f"Genomic position",
        fontsize=14
        )
    
    plt.ylabel(
        "B-allele frequency",
        fontsize=14
        )
    
    plt.title(
        f"B-allele frequency — chromosome {chromosome}",
        fontsize=16
        )

    # Get current axes
    ax = plt.gca()

    # -------------------------
    # X-axis
    # -------------------------

    # Major ticks every 10 Mb
    ax.xaxis.set_major_locator(
        MultipleLocator(1e7)
        )
    
    # Minor ticks every 1 Mb
    ax.xaxis.set_minor_locator(
        MultipleLocator(1e6)
        )
    
    # Format major tick labels as Mb
    ax.xaxis.set_major_formatter(
        FuncFormatter(
            lambda x, 
            pos: f"{x/1e6:.0f} Mb"
            )
            )
    # Make major and minor ticks visually distinct
    ax.tick_params(
        axis="x",
        which="major",
        length=7
    )

    ax.tick_params(
        axis="x",
        which="minor",
        length=4
    )

    # -------------------------
    # Y-axis
    # -------------------------
    # Major ticks every 0.2
    ax.yaxis.set_major_locator(
        MultipleLocator(0.2)
    )

    # Minor ticks every 0.1
    ax.yaxis.set_minor_locator(
        MultipleLocator(0.1)
    )

    # Y-axis tick sizes
    ax.tick_params(
        axis="y",
        which="major",
        length=7
    )

    ax.tick_params(
        axis="y",
        which="minor",
        length=4
    )

    # Grid
    plt.grid(
        True,
        alpha=0.3)
    
    plt.legend()

    # Save the plot
    plt.savefig(
        f"{sample_id}_chrom_{chromosome}.png", 
        dpi=300, 
        bbox_inches="tight"
        )
    # Close the figure before creating the next one
    plt.close()