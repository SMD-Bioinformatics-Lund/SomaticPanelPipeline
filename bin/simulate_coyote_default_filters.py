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

def main(args):
    """
    """
    config = _read_config_json(args.config)

    if args.known:
        known,known_statistics = _read_known_variants(args.known)
    else:
        known = {}

    statistics = _read_and_filter_vcf(args.vcf, config, args.tumor_id, known)

    if args.known:
        _print_summary(statistics, known_statistics)
    else:
        _print_summary(statistics)

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
    variants_missed = 0
    false_postives_found = 0
    tiered_found = 0

    known_summary = {}

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
        # for known variants
        is_match = False
        is_false = False
        
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


        # known
        if simple_id in known:
            is_match = _match_known(all_hgvsp,all_hgvsc,known[simple_id])
            if is_match:
                known_summary[simple_id] = "found"
                if is_match.get("fp") == "False positive" or is_match.get("irr") == "Irrelevant":
                    is_false = True
                    false_postives_found +=1
                elif is_match.get('Tier') <= TIER:
                    tiered_found +=1

        if var_is_kept and var_is_shown:
            variants_found +=1
        # any variants missed from known due to filters
        else:
            if is_match and not is_false:
                known_summary[simple_id] = reason_for_filter
    
    # if known was provided
    if len(known.keys()) > 0:
        for kvar in known:
            if kvar not in known_summary:
                for variant in known[kvar]:
                    if variant.get('Tier') <= TIER and (variant.get("fp") != "False positive" or variant.get("irr") != "Irrelevant"):
                        variants_missed +=1
                        known_summary[kvar] = "not called/tiered"
                    else:    
                        known_summary[kvar] = "not called/fp/irr"
        statistics['false_postives_found'] = false_postives_found
        statistics['tiered_variants_found'] = tiered_found
        statistics['tiered_variants_missed'] = variants_missed
    
    pprint(known_summary)
    statistics['variants_found_under_filters'] = variants_found
        
    return statistics

def _truthy(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}
    

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
    known_statistics = {
        'false_positives' : 0,
        'found' : 0,
        'tiered' : 0,
        'germline' : 0,
    }
    with open(filepath) as f:
        variants_ = csv.DictReader(f)
        for row in variants_:
            coordinates = _known_coordinates(row)

            # make Tier integers
            tier = row.get('Tier') or row.get('tier')
            try:
                tier = int(tier)
                if tier <= TIER:
                    known_statistics['tiered'] +=1
            except:
                tier = 10
            row['Tier'] = tier

            # save germline counts
            row['Flags'] = row.get('Flags') or row.get('flags') or ""
            if "GERM" in row.get('Flags'):
                known_statistics['germline'] +=1

            # make HGVS nomeclature more code friendly
            hgvs = row.get('HGVS') or ""
            hgvs_p = None
            hgvs_c = None
            p_match = re.search(r"(p\.[^| ]+)", hgvs)
            c_match = re.search(r"(c\.[^| ]+)", hgvs)
            if p_match:
                hgvs_p = p_match.group(1)
                row['hgvsp'] = hgvs_p
            elif row.get('hgvsp'):
                row['hgvsp'] = row['hgvsp'].strip()
            if c_match:
                hgvs_c = c_match.group(1)
                row['hgvsc'] = hgvs_c
            elif row.get('hgvsc'):
                row['hgvsc'] = row['hgvsc'].strip()

            if _truthy(row.get("false_positive")):
                row["fp"] = "False positive"
            if _truthy(row.get("irrelevant")):
                row["irr"] = "Irrelevant"

            if row.get("fp") == "False positive" or row.get("irr") == "Irrelevant":
                known_statistics['false_positives'] +=1
            
            variants[coordinates].append(row)

    known_statistics["found"] = len(variants.keys())
    return variants,known_statistics

def _known_coordinates(row):
    if row.get('Chr:Pos'):
        return row['Chr:Pos'].replace("'","").replace(":","_")
    if row.get('chr_pos'):
        return row['chr_pos'].replace("'","").replace(":","_")
    if row.get('chrom') and row.get('pos'):
        return f"{row.get('chrom')}_{row.get('pos')}"
    raise ValueError("known variants CSV requires Chr:Pos, chr_pos, or chrom/pos columns")

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

def _soft_filters_coyote(
        config: dict,
        statistics: dict,
        tumor_af: float,
        pop_af: float,
        variant_depth: int,
        depth: int,
        normal_af: float,
        all_consequences: list,
        reason_for_filter: list
    ):
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


def _print_summary(statistics, known_statistics=None):
    print("\n" + "="*60)
    print("FILTER SUMMARY")
    print("="*60)

    print("\nHard filter failures:")
    print(f"  PON: {statistics['fail_pon']}")
    print(f"  Other: {statistics['fail_other_filter']}")
    print(f"  Population (hard): {statistics['fail_pop_load']}")

    print("\nSoft filter failures:")
    print(f"  Tumor VAF: {statistics['fail_tvaf']}")
    print(f"  Population (soft): {statistics['fail_pop_filter']}")
    print(f"  Control VAF: {statistics['fail_control_vaf']}")
    print(f"  Depth: {statistics['fail_dp']}")
    print(f"  Alt reads: {statistics['fail_vd']}")
    print(f"  VEP consequence: {statistics['fail_vep']}")

    print("\nRetained overrides:")
    print(f"  Germline retained: {statistics['retained_germline']}")

    print("\nFinal counts:")
    print(f"  Variants passing all filters: {statistics['variants_found_under_filters']}")

    if known_statistics:
        print("\n" + "="*60)
        print("KNOWN VARIANT PERFORMANCE")
        print("="*60)

        print(f"  Known variants (input): {known_statistics['found']}")
        print(f"  Tier ≤ {TIER}: {known_statistics['tiered']}")
        print(f"  False positives in known: {known_statistics['false_positives']}")
        print(f"  Germline in known: {known_statistics['germline']}")

        print("\nDetection:")
        print(f"  Tiered variants found: {statistics['tiered_variants_found']}")
        print(f"  Tiered variants missed: {statistics['tiered_variants_missed']}")
        print(f"  Known false positives found: {statistics['false_postives_found']}")

        # optional recall metric
        denom = statistics['tiered_variants_found'] + statistics['tiered_variants_missed']
        if denom > 0:
            recall = statistics['tiered_variants_found'] / denom
            print(f"\n  Recall (tier ≤ {TIER}): {recall:.3f}")

if __name__ == "__main__":
    cli()
