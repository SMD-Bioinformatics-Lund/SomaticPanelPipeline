#!/usr/bin/env python3

"""
Description:
This script extracts informative SNPs from a posterior file based on LOH (Loss of Heterozygosity) regions provided in an input file (*_purecn_loh_filtered.csv).
These extracted SNPs are the ones used by PureCN for making LOH calls. For those extracted SNPs, a plot is generated for each chromosomes,
for the corresponding observed allele ratio (MAF) (or also called BAF).
The script also generates a bed files for GENS import for these CN_LOH calls.

Inputs:
1. --input_loh: A TSV file containing LOH calls with columns for chromosome, start, end, and LOH type.
2. --posterior_file: A TSV file containing SNP information with various attributes.

Outputs:
--out_snp: A CSV file containing SNPs that meet the filtering criteria and fall within the specified CN_LOH regions.
--out_gens (-o2) A TSV file for GENS visualization of the Copy Neutral LOH events (CN-LOH).
- Optional plot for each chromosome showing the B-allele frequency (BAF) for the extracted SNPs if CN_LOH calls are present (_baf_purecn_cnlohseg.png)
--out_coyote (-o3) A JSON file for Coyote visualization of the Copy Neutral LOH events (CN-LOH).

Usage Example:
cnloh_visualization_and_extractSNP.py \
    --input_loh sampleID_purecn_loh_filtered.csv \
    --posterior_file sampleID_snv_posteriors.tsv \
    --out_snp sampleID_informative_cnloh_snp.tsv \
    --out_gens sampleID_purecn_cnloh_gens.tsv \
    --out_coyote sampleID_cnloh_coyote.json
"""

import argparse
import csv
import json
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
        help='<sample>_purecn_loh_filtered.csv is used to collect the chromosome and positions of the LOH calls',
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
        help='tsv file _informative_cnloh_snp.tsv containing the extracted informative SNPs used by PureCN for making LOH calls',
        required=False
    )

    parser.add_argument(
        '-o2',
        '--out_gens',
        help='tsv file _purecn_cnloh_gens.tsv for visualization of the Copy Neutral LOH events (CN-LOH)',
        required=True
    )

    parser.add_argument(
        '-o3',
        '--out_coyote',
        help='json file _cnloh_coyote.json for visualization of the Copy Neutral LOH events (CN-LOH) in Coyote',
        required=False
    )

    return parser.parse_args()

def main():
    args = parse_args()

    # dictionary to collect LOH calls
    cnloh_regions = {} 
# ------------------------------------------------------------
# Generate the gens_tsv file with the CN-LOH segment
# for GENS visualization (violet segment)
# -----------------------------------------------------------
    coyote_json = []

    with open(args.input_loh, "r", newline="") as loh_call, \
        open(args.out_gens, "w", newline="") as gens_tsv:

        reader = csv.DictReader(loh_call)
        # write GENS header
        gens_tsv.write(
            "chromosome\tstart\tend\tname\tcomments\tcolor\n"
            )
        
        for row in reader:
            # Read values from the PureCN CSV
            chrom = row["chr"].strip()
            start = float(row["start"])
            end = float(row["end"])
            loh_type = row["type"].strip()

            snp = int(row["num.snps"])
            delta_maf = float(row["delta_maf"])
            snp_per_mb = float(row["snp_per_mb"])

            # Purity is numeric 
            purity = float(row["purity_purecn"])

            # Purity warning is either "." or "unreliable"
            purity_warning = row["purity_warning"].strip()

            # Keep information for CN-LOH calls
            if "COPY-NEUTRAL LOH" in loh_type:
                cnloh_regions.setdefault(chrom, []).append((start, end))

            # Information displayed by GENS
                comments = (
                    f"SNP: {snp};"
                    f"BAFexpected - BAF observed: {delta_maf};"
                    f"SNPs per Mb: {snp_per_mb};"
                    f"Purity: {purity};"
                    f"Purity accuracy warning: {purity_warning}"
                )

                # Write GENS TSV
                gens_tsv.write(
                    f"{chrom}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"CN-LOH\t"
                    f"{comments}\t"
                    f"rgb(127,0,255)\n"
                    )
                
                # Add entry for the Coyote Json
                coyote_json.append({
                "chromosome": int(chrom) if chrom.isdigit() else chrom,
                "start": int(start),
                "end": int(end),
                "name": "CN-LOH"
            })
    # Write JSON file
    with open(args.out_coyote, "w") as json_file:
        json.dump(coyote_json, json_file, indent=2)

    # -----------------------------------------------------------
    # Generate the file containing all informative SNPs for the CN-LOH calls used by PureCN
    # -----------------------------------------------------------
    with open(args.posterior_file, "r") as posterior, \
        open(args.out_snp, "w") as snps:

        for line in posterior:

            # Keep header
            if line.startswith("chr"):
                snps.write(line)
                continue

            field = line.split("\t")
            
            chrom = field[0].strip()
            position = int(field[1].strip())

            # Skip chromosomes without CN-LOH regions
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
            
            # -----------------------------------------------------------
            # Check whether SNP position is inside an LOH interval
            # -----------------------------------------------------------
            position_in_loh = any(
                start <= position <= end
                for start, end in cnloh_regions[chrom]
                )

            if not position_in_loh:
                continue

            # SNP passes all filters: full original row is kept for that SNP
            snps.write(line)

    # -----------------------------------------------------------
    # MAKING THE PLOTs for BAF for the final CN LOH calls (one plot per chromosome)
    # -----------------------------------------------------------
    informative_snps = pd.read_csv(args.out_snp, sep="\t")

    # Get sample ID from output filename to use it for the plots
    sample_id = os.path.basename(args.out_snp).replace("_informative_cnloh_snp.tsv", "")

    # Make one plot for each chromosome
    for chromosome, data in informative_snps.groupby("chr"):

        plt.figure(figsize=(14, 6))

        # plot BAF
        plt.scatter(
            data["start"], 
            data["AR"], 
            facecolors="hotpink", 
            label="BAF"
            )
        
        # plot BAF adjusted for mapping biais
        plt.scatter(
            data["start"], 
            data["AR.ADJUSTED"], 
            facecolors="darkgreen", 
            label="BAF(mapping-biais adjusted)"
            )

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
            f"{sample_id}_chr_{chromosome}_{int(start)}_{int(end)}_CNLOH_BAF.png", 
            dpi=300, 
            bbox_inches="tight"
            )
        # Close the figure before creating the next one
        plt.close()

if __name__ == "__main__": 
    main()