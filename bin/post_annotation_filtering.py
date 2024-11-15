from pysam import VariantFile
from pprint import pprint
import cmdvcf
import argparse
import fnmatch


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
    

def read_vcf(infile,filters,af_cutoff):
    """
    
    """
    vcf_object = VariantFile(infile)
    
    print(vcf_object.header,end="")
    
    for var in vcf_object.fetch():
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        
        af_pass = check_gnomad_vaf(var_dict,af_cutoff)
        
        filter_pass = filter_flags(var_dict,filters)
        
        if af_pass and filter_pass:
            print(str(var),end="")


def main():
    af_cutoff = 0.05
    filters = []
    
    parser = argparse.ArgumentParser(description="This is a script to filter variants before loading into coyote. This can be based upon gnomad frequencies and filter flags")

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
    args = parser.parse_args()
    if args.filters is not None:
        filters = args.filters.split(',')
    if args.max_freq is not None:
        af_cutoff = args.max_freq
    
    read_vcf(args.vcf,filters,af_cutoff)
    

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