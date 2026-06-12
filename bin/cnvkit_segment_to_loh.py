#!/usr/bin/env python3

"""
Description: This script identifies LOH events and classifies their type based on the input cnvkitsegment (cnvkit_hrd) 
from CNVKIT_CALL (or CNVKIT_CALL_TC). 
This python script outputs these LOH events with category labels (tab delimited file) as well as a corresponding bed file for visualization.

Procedure:  1) filters on the minor allele (cn2) == 0 (LOH)
            2) Attributes a LOH type based on total_cn (total allele copy number)
            3) Flags probable LOH when both cn1 and cn2 are missing but total_cn == 1 
            (single copy loss with no BAF data is suspicious for LOH but cannot be confirmed
"""

import argparse

parser = argparse.ArgumentParser(description="Process CNVkit .cns output to identify and classify LOH events.")
parser.add_argument('--input',   required=True, help='input .cns file from CNVkit (tab-delimited) in CNVKIT_CALL(_TC) process')
parser.add_argument('--sample',  required=True, help='sample ID to use in the output loh_cat')
parser.add_argument('--sex',     required=True, help='sex to drop Y chr detection if female')
parser.add_argument('--out_cat', required=True, help='output filename for LOH results')

args = parser.parse_args()

with open(args.input, "r") as cns, \
     open(args.out_cat, 'w') as out_file:

    out_file.write("Chromosome\tStart_position\tEnd_position\ttotal_cn\tmajor_cn\tminor_cn\tloh_type\tsampleID\n")

    header = next(cns)
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
        cn1_raw = fields[9]   # major allele cn1, cn1 >= cn2 guaranteed by CNVkit > v0.9.0
        cn2_raw = fields[10]  # minor allele cn2

        # Skip Y chromosome if meta indicates a female patient
        if args.sex == "F" and chromosome in {"chrY", "Y"}:
            continue

        # Parse total CN, required, skip if absent or non-integer
        try:
            total_cn = int(cn_raw)
        except ValueError:
            continue
        
        # Parse allele-specific CNs
        # if either cn1 or cn2 is missing, we can deduct the missing from total_cn.
        # if cn1 and cn2 is missing, we cannot confirm LOH, but we can still flag probable LOH if total_cn=1
        cn1_empty = cn1_raw.strip() == ""
        cn2_empty = cn2_raw.strip() == ""

        if not cn1_empty and not cn2_empty:
            # Both present — normal case
            try:
                major_cn = int(cn1_raw)
                minor_cn = int(cn2_raw)
            except ValueError:
                continue
            both_missing = False

        elif cn1_empty and cn2_empty:
            # Both missing — cannot determine allele state
            major_cn = None
            minor_cn = None
            both_missing = True

        elif cn1_empty:
            # cn1 (major_cn) missing only, deduct it from total
            try:
                minor_cn = int(cn2_raw)
                major_cn = total_cn - minor_cn
            except ValueError:
                continue
            both_missing = False

        else:
            # cn2 (minor_cn) missing only, deduct it from total
            try:
                major_cn = int(cn1_raw)
                minor_cn = total_cn - major_cn
            except ValueError:
                continue
            both_missing = False


        # LOH classification
        if both_missing is False and minor_cn == 0:
            # Confirmed LOH: minor allele is absent
            if total_cn == 0:
                loh_type = "HOM_DEL" # homozygous deletion, both alleles lost
            elif total_cn == 1:
                loh_type = "LOH_DEL" # LOH with single copy loss, one allele lost
            elif total_cn == 2:
                loh_type = "CN_LOH" # copy-neutral LOH, one allele lost and the other duplicated (total CN=2 but minor CN=0)
            elif total_cn > 2:
                loh_type = "LOH_AMP" # LOH with amplification, one allele lost but total CN > 2 (e.g. one allele lost but the other is amplified)
            else:
                loh_type = "undetermined"

        elif both_missing and total_cn == 1:
            # Probable LOH: single copy but allele-specific CN unavailable
            # total_cn=1 with no BAF data cannot be confirmed, but is suspicious
            loh_type = "PROB_LOH_DEL"

        else:
            continue  # not LOH, skip

        # Format None as '.' for TSV output
        major_str = str(major_cn) if major_cn is not None else "."
        minor_str = str(minor_cn) if minor_cn is not None else "."

        out_file.write(f"{chromosome}\t{start}\t{end}\t{total_cn}\t{major_str}\t{minor_str}\t{loh_type}\t{args.sample}\n")
