#!/usr/bin/env python3

"""
Description: this script uses as inputs:
- the PureCN log file (to collect warning of reliability of the purity estimate)
- PureCN purity csv file (to collect purity estimate)
- PureCN loh_file (all LOH calls from PureCN)
 and discards irrelevant CNV calls.
Delta MAF (minor allele frequency) between the modelled and the observed value is calculated, as well as the number of SNPs per Mb. 
One output is generated, containing the filtered LOH calls with informative SNPs, M.flagged == FALSE and type with no empty or NA values, 
as well as the calculated delta MAF and SNPs per Mb, the purity estimate and a warning if the purity estimate is unreliable.

Usage:
filter_purecn_loh_output.py -l <log_file> -p <purity_file> -i <loh_file> -o <output_file>

"""
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Creating final LOH output table"
    )

    parser.add_argument(
        '-l',
        '--log_file',
        help='log file to collect potential warning of purity estimate unreliability',
        required=True
    )

    parser.add_argument(
        '-p',
        '--purity_file',
        help='csv file where final selected purity estimate is indicated',
        required=True
    )

    parser.add_argument(
        '-i',
        '--loh_file',
        help='loh csv with all LOH calls to be filtered',
        required=True
    )

    parser.add_argument(
        '-o',
        '--output_file',
        help='_loh_filtered.tsv containing LOH calls with informative SNPs, M.flagged == FALSE and type with no empty or NA values',
        required=True
    )

    return parser.parse_args()

def main():

    args = parse_args()
    
    # Check PureCN warning in the log file
    warning = "."
    with open(args.log_file, "r") as log_input:
        for line in log_input:
            if "No copy number alterations found, purity estimate unreliable." in line:
                print(f"From PureCN log: {line.strip()}")
                warning = "unreliable"
                break    
    
    # get estimated purity by PureCN in the csv sample output file
    purity = "."
    with open(args.purity_file, "r") as purity_input:
        lines = purity_input.readlines()
        for line in lines:
            if line.startswith("Sampleid"):
                continue
            if not line.startswith("Sampleid"):
                field = line.split(",")
                purity = field[1].strip()

    # Filter LOH calls from the sample_loh.csv output file
    with open(args.loh_file, "r") as loh_input, open(args.output_file, "w") as loh_filt:
        lines = loh_input.readlines()
        for line in lines:
            if line.startswith("Sampleid"):
                loh_filt.write(line.strip() + ",delta_maf,snp_per_mb,purity_purecn,purity_warning\n")
            if not line.startswith("Sampleid"):
                field = line.split(",")
                start = field[2].strip()
                end = field[3].strip()
                loh_type = field[7].strip()
                snps = field[10].strip()
                minor_flag = field[11].strip()
                maf_exp = field[12].strip()
                maf_obs = field[13].strip()

                if maf_obs == "NA":
                    continue

                # Calculate SNPs per Mb
                snp_per_mb = round(float(snps) / ((float(end) - float(start)) / 1000000), 2)

                # Calculate delta MAF
                if maf_exp == "NA":
                    delta_maf = "NA"

                else:
                    delta_maf = round(abs(float(maf_exp) - float(maf_obs)), 3)
                    
                # Keep relevant LOH calls
                if loh_type != "" and loh_type != "NA" and minor_flag == "FALSE" and int(snps) >= 10:
                        loh_filt.write(line.strip() + f",{delta_maf},{snp_per_mb},{purity},{warning}\n")

if __name__ == "__main__":
    main()
