#!/usr/local/bin/python
from pysam import VariantFile
from pprint import pprint
import cmdvcf
import argparse
import fnmatch
import re


def filter_flags(var,filters):
    """
    Pattern match exclusion filters, works with wildcards. If matched_items above 0 dont print variant
    """
    matched_items = 0
    var_filters = var["FILTER"].split(";")
    for pattern in filters:
        matched_items = matched_items + len([item for item in var_filters if fnmatch.fnmatch(item, pattern)])
    
    if matched_items > 0:
        return False
    else:
        return True
    

def override_filters(var, override_filters) -> bool:
    """
    If a variant is tagged as failed but keeps it if additional key, value terms are provided for exmaple is `likely pathogenic and pathogenic` in CLIN_SIG are kept
    """
    for filter_key_values in override_filters:
        for filter_key, filter_values in filter_key_values.items():
            if filter_key in var["INFO"]["CSQ"][0]:
                if any([value in var["INFO"]["CSQ"][0][filter_key] for value in filter_values]):
                    return True, filter_key
    return False, None

def read_vcf(infile,filters,af_cutoff, override_key_values):
    """
    
    """
    vcf_object = VariantFile(infile)
    
    print(vcf_object.header,end="")
    
    for var in vcf_object.fetch():
        override = False
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        
        af_pass = check_gnomad_vaf(var_dict,af_cutoff)
        
        filter_pass = filter_flags(var_dict,filters)

        if not filter_pass:
            override, override_key = override_filters(var_dict, override_key_values)

        if override:
            current_filter_column = var_dict["FILTER"]

            if "PASS" not in current_filter_column.split(";"):
                new_filter_column = f"PASS;{override_key}_Override;{current_filter_column}"
            else:
                new_filter_column = f"{override_key}_Override;{current_filter_column}"

            var = re.sub(current_filter_column, new_filter_column, str(var))
        
        if af_pass and (filter_pass or override):
            print(str(var),end="")


def main():
    af_cutoff = 0.05
    filters = []
    override_key_values = []
    
    parser = argparse.ArgumentParser(description="his is a script to filter variants before loading into coyote. This can be based upon gnomad frequencies and filter flags. This will not filter out variants that are tagged as failed but are `likely pathogenic and pathogenic`.")

    parser.add_argument(
        "--vcf",
        "-f",
        type=str,
        required=True,  # Makes this flag mandatory
        help="Path to the file to process"
    )
    parser.add_argument(
        "--max_freq",
        type=float,
        help="gnomAD frequency cutoff default 0.05"
    )
    parser.add_argument(
        "--filters",
        type=str,
        help="comma-separated list of filters to NOT keep"
    )
    parser.add_argument(
        "--override_key_values",
        type=str,
        help="semiclon-separated list of terms to keep if the variant is not passed by the filters, eg: 'CLIN_SIG=likely_pathogenic,pathogenic;CANONICAL=Yes,Y'"
    )
    args = parser.parse_args()
    if args.filters is not None:
        filters = args.filters.split(',')
    if args.max_freq is not None:
        af_cutoff = args.max_freq

    if args.override_key_values:
        override_key_values = [
            {kv.split('=')[0]: kv.split('=')[1].split(',')}
            for kv in args.override_key_values.split(';')
            if kv
        ]
    
    read_vcf(args.vcf,filters,af_cutoff, override_key_values)
    

def check_gnomad_vaf(var,af_cutoff):
    """
    Find gnomad annotations, return true or false if freq is above cutoff
    """
    af_dict           = {}
    gnomad            = var["INFO"]["CSQ"][0].get("gnomAD_AF", 0)
    gnomad_genome     = var["INFO"]["CSQ"][0].get("gnomADg_AF", 0)
    gnomad_max        = var["INFO"]["CSQ"][0].get("MAX_AF", 0)

    af = 0
    if gnomad:
        gnomad = max_gnomad(gnomad)
        af = gnomad
        if gnomad == '':
            af = 0
    elif gnomad_genome:
        gnomad = max_gnomad(gnomad)
        af = gnomad
        if gnomad == '':
            af = 0
    else:
        ## could not find gnomad annotations
        return 1
    
    if float(af) >= af_cutoff:
        return 0
    else:
        return 1

def max_gnomad(gnomad):
    """
    check if gnoamd is multivalued, split and max
    """
    try:
        gnomad_list = gnomad.split('&')
        if gnomad_list:
            return float(max(gnomad_list))
    except:
        return gnomad

if __name__ == "__main__":
    main()