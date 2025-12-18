#!/usr/bin/env python
import logging
import argparse
import json
import cmdvcf
import csv
import re

from pysam import VariantFile
from collections import defaultdict
from pprint import pprint

GNOMAD_KEY = "gnomAD_AF"
TIER = 3
#known_variants_missed.append(f"Pop_af:{pop_af}  tumor_af:{tumor_af}  variant_depth:{variant_depth} depth:{depth} normal_af:{normal_af}  filters:{filters} vep_consequences:{all_consequences}")

def main(args):
    """
    """
    config = _read_config_json(args.config)

    if args.known:
        known = _read_known_variants(args.known)
    else:
        known = {}

    statistics = _read_and_filter_vcf(args.vcf,config,args.tumor_id,known)
    pprint(statistics)

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
    variants_found = 0
    false_postives_found = 0
    tiered_found = 0
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
        
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        simple_id = f"{var.chrom}_{var.pos}"

        filters = var_dict['FILTER']
        pop_af = max_gnomad(var_dict['INFO']['CSQ'][0][GNOMAD_KEY])
        if pop_af == "":
            pop_af = 0
        tumor_af, variant_depth, depth, normal_af = _get_af(var_dict,tumor_id)
        all_consequences,all_hgvsp,all_hgvsc = _get_all_consequences(var_dict['INFO']['CSQ'])

        reason_for_filter = []
        var_is_kept = _hard_filters_cli(config,filters,statistics,pop_af,reason_for_filter)
        var_is_shown = _soft_filters_coyote(config,statistics,tumor_af,pop_af,variant_depth, depth, normal_af,all_consequences,reason_for_filter)

        # completely override all rules for this mf
        if "GERMLINE" in filters and not 'GERMLINE_RISK' in filters:
            var_is_kept = True
            var_is_shown = True
            statistics['retained_germline'] +=1

        if var_is_kept and var_is_shown:
            variants_found +=1
            if simple_id in known:
                is_match = _match_known(all_hgvsp,all_hgvsc,known[simple_id])
                if is_match:
                    if is_match.get("fp") == "False positive" or is_match.get("irr") == "Irrelevant":
                        false_postives_found +=1
                        continue
                    if is_match.get('Tier') <= TIER:
                        tiered_found +=1
        else:
            if simple_id in known:
                ...
        
    statistics['variants_found_under_filters'] = variants_found
    statistics['false_postives_found'] = false_postives_found
    statistics['tiered_found'] = tiered_found
    return statistics
    

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
    all_hgvsp = []
    all_hgvsc = []
    for trans in csq:
        for quence in trans['Consequence']:
            all_consequences.append(quence)
        hgvsp = trans.get("HGVSp")
        if isinstance(hgvsp, str):
            all_hgvsp.append(hgvsp.split(":", 1)[-1].strip())
        hgvsc = trans.get("HGVSc")
        if isinstance(hgvsc, str):
            all_hgvsc.append(hgvsc.split(":", 1)[-1].strip())
    all_consequences = list(set(all_consequences))
    return all_consequences,all_hgvsp,all_hgvsc

def _read_known_variants(filepath):
    variants = defaultdict(list)
    with open(filepath) as f:
        variants_ = csv.DictReader(f)
        for row in variants_:
            coordinates = row['Chr:Pos'].replace("'","").replace(":","_")
            tier = row.get('Tier')
            try:
                tier = int(tier)
            except:
                tier = 10
            row['Tier'] = tier
            hgvs = row['HGVS']
            hgvs_p = None
            hgvs_c = None
            p_match = re.search(r"(p\.[^| ]+)", hgvs)
            c_match = re.search(r"(c\.[^| ]+)", hgvs)
            if p_match:
                hgvs_p = p_match.group(1)
                row['hgvsp'] = hgvs_p
            if c_match:
                hgvs_c = c_match.group(1)
                row['hgvsc'] = hgvs_c
            
            variants[coordinates].append(row)

    return variants

def _hard_filters_cli(config,filters,statistics,pop_af,reason_for_filter):
    """
    """
    var_is_kept = True
    # hard filters, dont make it past CLI
    for filter_name in config['cli_filters']['filters']:
        if filter_name in filters:
            if 'PON' in filter_name:
                statistics['fail_pon'] +=1
                reason_for_filter.append("PON")
            else:
                statistics['fail_other_filter'] +=1
                reason_for_filter.append("failed_other_filter")
            var_is_kept = False
    if pop_af > config['cli_filters']['af']:
        var_is_kept = False
        statistics['fail_pop_load'] +=1
        reason_for_filter.append("hard_popfreq")

    return var_is_kept

def _soft_filters_coyote(config,statistics,tumor_af,pop_af,variant_depth, depth, normal_af,all_consequences,reason_for_filter):
    """
    """
    var_is_shown = True
    # soft filters, what is shown by default in coyote
    if tumor_af >= config['snv_filters']['max_freq'] or tumor_af <= config['snv_filters']['min_freq']:
        var_is_shown = False
        statistics['fail_tvaf'] +=1
        reason_for_filter.append("tvaf")
    if pop_af > config['snv_filters']['max_popfreq']:
        var_is_shown = False
        statistics['fail_pop_filter'] +=1
        reason_for_filter.append("soft_popfrq")
    if depth < config['snv_filters']['min_depth']:
        var_is_shown = False
        statistics['fail_dp'] +=1
        reason_for_filter.append("min_depth")
    if variant_depth < config['snv_filters']['min_alt_reads']:
        var_is_shown = False
        statistics['fail_vd'] +=1
        reason_for_filter.append("min_alt_reads")
    if normal_af is not None:
        if normal_af >= config['snv_filters']['max_control_freq']:
            var_is_shown = False
            statistics['fail_control_vaf'] +=1
            reason_for_filter.append("nvaf")
    if not (set(all_consequences) & set(config['vep_filters'])):
        var_is_shown = False
        statistics['fail_vep'] +=1
        reason_for_filter.append("vep_consequence")

    return var_is_shown

def _match_known(all_hgvsp,all_hgvsc,matches):
    """
    """
    for match in matches:
        if match.get('hgvsp') in all_hgvsp:
            return match
        if match.get('hgvsc') in all_hgvsc:
            return match
    return False

if __name__ == "__main__":
    cli()