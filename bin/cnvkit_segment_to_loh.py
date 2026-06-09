#!/usr/bin/env python3

"""
Description: This script identifies LOH events and classifies their type based on the output cnvkitsegment from CNVKIT_CALL (or CNVKIT_CALL_TC). 
It outputs these LOH events as well as a corresponding bed file for visualization in Coyote.

Necessary input is the txt output file from scarHRD.

Procedure:  1) filter on the minor allele (cn2) == 0 (LOH)
            2) Attributes a LOH type based on total_cn (total allele copy number)
"""

import argparse

parser = argparse.ArgumentParser(description="Process CNVkit .cns output to identify and classify LOH events.")
parser.add_argument('--input',   required=True, help='input .cns file from CNVkit (tab-delimited)')
parser.add_argument('--sample',  required=True, help='sample ID to use in output')
parser.add_argument('--out',     required=True, help='output filename for LOH results')
parser.add_argument('--out_bed', required=True, help='output filename for LOH BED file')
args = parser.parse_args()

with open(args.input, "r") as cns, \
     open(args.out,     'w') as out_file, \
     open(args.out_bed, 'w') as bed_file:

    out_file.write("sampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tmajor_cn\tminor_cn\tloh_type\n")
    bed_file.write("Chromosome\tStart_position\tEnd_position\tsampleID\tloh_type\n")

    header = next(cns, None)
    if not header:
        raise ValueError("Input file is empty or missing a header.")

    for line in cns:
        fields = line.strip().split("\t")

        chromosome = fields[0]
        start      = fields[1]
        end        = fields[2]
        # fields[3] = gene, fields[4] = log2, fields[5] = baf (may be empty)
        # fields[6] = ci_hi, fields[7] = ci_lo
        cn_raw  = fields[8]
        cn1_raw = fields[9]   # major allele, cn1 >= cn2 guaranteed by CNVkit > v0.9.0
        cn2_raw = fields[10]  # minor allele

        # Skip segments where CN calls are missing or non-integer
        try:
            total_cn  = int(cn_raw)
            major_cn  = int(cn1_raw)
            minor_cn  = int(cn2_raw)
        except ValueError:
            continue

        # LOH defined as minor allele == 0
        if minor_cn == 0:
            if total_cn == 0:
                loh_type = "HOM_DEL"
            elif total_cn == 1:
                loh_type = "LOH_DEL"
            elif total_cn == 2:
                loh_type = "CN_LOH"
            elif total_cn > 2:
                loh_type = "LOH_AMP"
            else:
                loh_type = "undetermined"

            out_file.write(f"{args.sample}\t{chromosome}\t{start}\t{end}\t{total_cn}\t{major_cn}\t{minor_cn}\t{loh_type}\n")
            bed_file.write(f"{chromosome}\t{start}\t{end}\t{args.sample}\t{loh_type}\n")