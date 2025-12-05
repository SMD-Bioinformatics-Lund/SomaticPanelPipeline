#!/usr/bin/env python
import logging
import argparse
import json
import cmdvcf

from pysam import VariantFile
from collections import defaultdict
from pprint import pprint

GNOMAD_KEY = "gnomAD_AF"

def main(args):
    """
    """
    config = _read_config_json(args.config)
    known = _read_known_variants(args.known)
    _read_and_filter_vcf(args.vcf,config,args.tumor_id,known)


def cli():
    """
    """
    parser = argparse.ArgumentParser(description="Evaluate number of shown SNV/indel variants under provided coyote default settings")
    parser.add_argument("--config", type=str, help="Path to assay JSON config")
    parser.add_argument("--vcf", type=str, help="Path to vep-annotated vcf")
    parser.add_argument("--known", type=str, help="Path to known variants")
    parser.add_argument("--tumor_id", type=str, help="sample id of tumor")
    args = parser.parse_args()
    if not args.config:
        exit("no config provided")
    main(args)


def _read_config_json(json_path):
    with open(json_path) as jf:
        configs = json.load(jf)
    return configs

def _read_and_filter_vcf(vcf_path,config,tumor_id,known):
    vcf_object = VariantFile(vcf_path)
    original_header = vcf_object.header
    variants_found_under_filter_settings = 0
    known_variants_missed = []

    statistics = {
        'fail_pon' : 0,
        'fail_other_filter' : 0,
        'fail_pop_load' : 0,
        'fail_control_vaf' : 0,
        'fail_pop_filter' : 0,
        'fail_tvaf' : 0,
        'fail_dp' : 0,
        'fail_vd' : 0,
        'retained_germline' : 0,
        'fail_vep' : 0
    }

    for var in vcf_object.fetch():
        var_is_kept = True
        var_is_shown = True
        should_be_found = False

        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        simple_id = f"{var.chrom}_{var.pos}_{var.ref}_{var.alts[0]}"
        even_simpler_id = f"{var.chrom}:{var.pos}"

        if even_simpler_id in known:
            should_be_found = True

        filters = var_dict['FILTER']
        pop_af = max_gnomad(var_dict['INFO']['CSQ'][0][GNOMAD_KEY])
        if pop_af == "":
            pop_af = 0
        tumor_af, variant_depth, depth, normal_af = _get_af(var_dict,tumor_id)
        all_consequences = _get_all_consequences(var_dict['INFO']['CSQ'])

        # hard filters, dont make it past CLI
        for filter_name in config['cli_filters']['filters']:
            if filter_name in filters:
                if 'PON' in filter_name:
                    statistics['fail_pon'] +=1
                else:
                    statistics['fail_other_filter'] +=1
                var_is_kept = False
        if pop_af > config['cli_filters']['af']:
            var_is_kept = False
            statistics['fail_pop_load'] +=1
        
        # soft filters, what is shown by default in coyote
        if tumor_af >= config['snv_filters']['max_freq'] or tumor_af <= config['snv_filters']['min_freq']:
            var_is_shown = False
            statistics['fail_tvaf'] +=1
        if pop_af > config['snv_filters']['max_popfreq']:
            var_is_shown = False
            statistics['fail_pop_filter'] +=1
        if depth < config['snv_filters']['min_depth']:
            var_is_shown = False
            statistics['fail_dp'] +=1
        if variant_depth < config['snv_filters']['min_alt_reads']:
            var_is_shown = False
            statistics['fail_vd'] +=1
        if normal_af is not None:
            if normal_af >= config['snv_filters']['max_control_freq']:
                var_is_shown = False
                statistics['fail_control_vaf'] +=1
        if not (set(all_consequences) & set(config['vep_filters'])):
            var_is_shown = False
            statistics['fail_vep'] +=1

        # completely override all rules for this mf
        if "GERMLINE" in filters and not 'GERMLINE_RISK' in filters:
            var_is_kept = True
            var_is_shown = True
            statistics['retained_germline'] +=1

        if var_is_kept and var_is_shown:
            variants_found_under_filter_settings +=1
            # else:
            #     print(f"{simple_id} should not be found")
            #     print(f"Pop_af:{pop_af}  tumor_af:{tumor_af}  variant_depth:{variant_depth} depth:{depth} normal_af:{normal_af}  filters:{filters} vep_consequences:{all_consequences}")
        else:
            if should_be_found:
                known_variants_missed.append(f"Pop_af:{pop_af}  tumor_af:{tumor_af}  variant_depth:{variant_depth} depth:{depth} normal_af:{normal_af}  filters:{filters} vep_consequences:{all_consequences}")

    print(f"variants shown in coyote:{variants_found_under_filter_settings} ")
    print(f"variants missed:{known_variants_missed} ")
    pprint(statistics)

def max_gnomad(gnomad):
    """
    check if gnoamd is multivalued, split and max
    """
    try:
        gnomad_list = gnomad.split("&")
        if gnomad_list:
            return float(max(gnomad_list))
    except:
        return gnomad

def _get_af(var_dict,tumor_id):
    tumor_af = None
    variant_depth = None
    depth = None
    normal_af = None
    for sample in var_dict['GT']:
        if sample['_sample_id'] == tumor_id:
            tumor_af = sample['VAF']
            variant_depth = sample['VD']
            depth = sample['DP']
        else:
            normal_af = sample['VAF']
    return tumor_af, variant_depth, depth, normal_af

def _get_all_consequences(csq):
    all_consequences = []
    for trans in csq:
        for quence in trans['Consequence']:
            all_consequences.append(quence)
    all_consequences = list(set(all_consequences))
    return all_consequences

def _read_known_variants(filepath):
    variants = []
    with open(filepath) as f:
        for line in f:
            variants.append(line.rstrip())
    return variants

if __name__ == "__main__":
    cli()