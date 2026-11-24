#!/usr/bin/env python3

"""
Description:
This script modifies the format of the VCF [output of MARK_GERMLINES process (SPP)], to be compatible with the input vcf for PureCN.
Modifications include:

- AD FORMAT field added: constructed as (DP-VD),VD for each sample column
  from the vcf:
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  ##FORMAT=<ID=VD,Number=1,Type=Integer,Description="ALT allele observation count">
  Organized as follow: GT:VAF:VD:DP

- POP_AF INFO field added: gnomADg_AF extracted from VEP CSQ

For POP_AF:
-gnomADg_AF is extracted from the CSQ field that looks at follow:
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. \
# Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc\
# |HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE\
# |STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS\
# |CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF\
# |BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF\
# |EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF\
# |gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS\
# |HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_PHRED|CADD_RAW|LoFtool|gnomADg\
# |gnomADg_AF_popmax|gnomADg_AF|gnomADg_popmax|COSMIC|COSMIC_CNT">

-If several transcript annotations have gnomADg_AF values, they must all be identical
-If conflicting values are found, the script stops with an error

Usage:
format_vcf_for_purecn_AD_POPAF.py --input <markgerm.vcf> --output <purecn.vcf>

"""
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add AD FORMAT and POP_AF INFO fields to VCF for PureCN"
    )

    parser.add_argument(
        '-i',
        '--input',
        help='Input vcf to be formatted (*.agg.pon.vep.markgerm.dbsnp.vcf)',
        required=True
    )

    parser.add_argument(
        '-o',
        '--output',
        help='formatted vcf for PureCN (*.purecn.vcf)',
        required=True
    )

    return parser.parse_args()


def construct_ad(format_string, sample_string):
    """
    Build an AD value (ref_count,alt_count) from DP and VD in a sample column.

    DP = total read depth
    VD = alt allele observation count
    ref_count = max(0, DP - VD)
    Returns '.,.' if DP or VD are missing or non-numeric.
    """

    format_keys = format_string.split(":")
    sample_vals = sample_string.split(":")
    format_map = dict(zip(format_keys, sample_vals))

    try:
        dp = int(format_map["DP"])
        vd = int(format_map["VD"])

        if dp < 0 or vd < 0:
            return ".,."

        ref_count = max(0, dp - vd)

        return f"{ref_count},{vd}"

    except (KeyError, ValueError):
        return ".,."


def extract_gnomad_af(csq_string, gnomad_af_idx):
    """
    Extract gnomADg_AF from the CSQ field.

    All non-missing gnomADg_AF values must be identical.
    Returns None if no gnomADg_AF value is present.
    """
    af_values = []  # collect all non-missing gnomADg_AF values from all transcripts for the variant
    
    for transcript in csq_string.split(","):  # a comma separates CSQ annotations for each transcript
        fields = transcript.split("|")  # all fields for one transcript annotation are separated by a pipe in the CSQ field

        if gnomad_af_idx < len(fields):
            value = fields[gnomad_af_idx]

            if value not in ("", "."):
                af_values.append(float(value))

    if not af_values:
        return None

    if len(set(af_values)) > 1:  # then different gnomADg_AF values were found for the same VCF variant
        raise ValueError(
            f"Conflicting gnomADg_AF values: {af_values}"
        )

    return af_values[0]


def main():
    """
    Main function to process the input VCF file, add AD FORMAT and POP_AF INFO fields,
    and write the formatted output VCF file for PureCN.
    """

    args = parse_args()

    with open(args.input, 'r') as in_vcf, open(args.output, 'w') as out_vcf:

        # Find the gnomADg_AF position in the CSQ header
        gnomad_af_idx = None

        for line in in_vcf:

            # Extracting the index of the gnomADg_AF field from the CSQ header line
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r'Format: ([^">]+)', line)  # r for raw string, thus pipes are not interpreted
                # one or more (+) character that is not " or > "(those latter are contained in the CSQ header line)

                if match:
                    csq_fields = match.group(1).split("|")  # group(1) returns the first match corresponding to [...], group(0) would return all the regex match

                    if "gnomADg_AF" in csq_fields:
                        gnomad_af_idx = csq_fields.index("gnomADg_AF")

                out_vcf.write(line)

                continue

            if line.startswith("#CHROM"):

                if gnomad_af_idx is None:
                    print(
                        "WARNING: gnomADg_AF not found in CSQ header. "
                        "POP_AF will not be added."
                    )
                else:
                    # write the new INFO header line for POP_AF
                    out_vcf.write(
                        '##INFO=<ID=POP_AF,Number=1,Type=Float,'
                        'Description="Population allele frequency from '
                        'gnomAD genomes (gnomADg_AF) extracted from VEP CSQ '
                        'for PureCN">\n'
                    )

                # write the new FORMAT header line for AD
                out_vcf.write(
                    '##FORMAT=<ID=AD,Number=R,Type=Integer, '
                    'Description="Allelic depths for ref and alt alleles (ref_count,alt_count), '
                    'derived from DP and VD for PureCN">\n'
                )

                out_vcf.write(line)

                # calculate how many samples are included in the vcf
                fields = line.strip().split('\t')
                num_samples = len(fields) - 9  # 9 columns are fixed, the rest are sample columns
                print(f'Number of samples detected in the vcf is {num_samples}')  # simple check of the vcf structure

                continue

            if line.startswith('#'):
                # Write other header lines unchanged
                out_vcf.write(line)

                continue

            # ---------------------------------------------------------
            # Process variant
            # ---------------------------------------------------------

            fields = line.strip().split('\t')

            chrom = fields[0]
            pos = fields[1]
            id_field = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_field = fields[6]
            info = fields[7]
            format_field = fields[8]
            sample_cols = fields[9:]  # 1 column tumor-only, 2 columns paired

            # ---------------------------------------------------------
            # Extract gnomADg_AF from CSQ in info
            # ---------------------------------------------------------

            pop_af = None

            if gnomad_af_idx is not None:
                csq_match = re.search(r"(?:^|;)CSQ=([^;]+)", info)
                # "(?:^|;)" => non-capturing group (?:) matching the start (^) of the string or (|) a semicolon (;), ensuring we capture the CSQ field correctly
                # the segment of interest starts with "CSQ=" and continues until the next semicolon or the end of the string
                # "([^;]+)" => capturing group that matches one or more characters that are not a semicolon: content of the CSQ field

                if csq_match:
                    try:
                        pop_af = extract_gnomad_af(
                            csq_match.group(1),
                            gnomad_af_idx
                        )

                    except ValueError as e:
                        print(
                            f"WARNING: {e}. "
                            f"Skipping variant at position {pos} "
                            f"on chromosome {chrom}."
                        )
                        pop_af = None

            # Add POP_AF to INFO
            if pop_af is not None:
                if info == ".":
                    info = f"POP_AF={pop_af:g}"
                else:
                    info += f";POP_AF={pop_af:g}"

            # ---------------------------------------------------------
            # Construct AD
            # ---------------------------------------------------------

            new_format = format_field + ":AD"

            new_sample_cols = [
                sample + ":" + construct_ad(format_field, sample)
                for sample in sample_cols
            ]

            out_fields = (
                [chrom, pos, id_field, ref, alt, qual, filter_field, info, new_format]
                + new_sample_cols
            )

            out_vcf.write("\t".join(out_fields) + "\n")


if __name__ == "__main__":
    main()