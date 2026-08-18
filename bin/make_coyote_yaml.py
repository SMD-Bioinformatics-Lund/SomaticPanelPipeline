#!/usr/bin/env python3
"""
Description: this script generates a YAML document for Coyote3 data import, combining meta info, optional file locations and other
parameters. This script was created to replace the currently used processes "COYOTE" and "COYOTE_YAML".
"""

import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate a COYOTE YAML document from SomaticPanelPipeline outputs and all_samples_metadata'
    )
    parser.add_argument('--group',            required=True, help='sample group name')
    parser.add_argument('--meta',             required=True, help='path to meta JSON file - one dict per sample')
    parser.add_argument('--vcf',              required=True, help='VCF filename')
    parser.add_argument('--json_info',        required=True, help='path to JSON file containing output filenames')
    parser.add_argument('--subdir',           required=True, help='pipeline subdirectory for building /access/ paths')
    parser.add_argument('--assay',            required=True, help='coyote group/assay name')
    parser.add_argument('--environment',      required=True, help='production, development, validation, or testing')
    parser.add_argument('--pipeline_name',    required=True, help='pipeline name from workflow manifest')
    parser.add_argument('--pipeline_version', required=True, help='pipeline version from workflow manifest')
    parser.add_argument('--out',              required=True, help='output YAML filename')
    return parser.parse_args()

def build_access_path(subdir, subfolder, filename):
    """Build the /access/ path used by Coyote to locate a file."""
    return f"/access/{subdir}/{subfolder}/{filename}"

def write_yaml(coyote_doc, out_filepath):
    """
    Write a flat YAML file.
    Handles None, bool, int/float and string values with correct YAML formatting.
    """
    with open(out_filepath, 'w') as yaml_file:
        yaml_file.write('---\n')  # Start of YAML document
        for yaml_key, yaml_value in coyote_doc.items():
            if yaml_value is None:
                 # YAML null for absent values
                yaml_file.write(f"{yaml_key}: null\n")
            elif isinstance(yaml_value, bool):
                # Checks if the value is a Python boolean
                # YAML booleans must be lowercase true/false
                yaml_file.write(f"{yaml_key}: {'true' if yaml_value else 'false'}\n")
            elif isinstance(yaml_value, (int, float)):
                # numbers are written without quotes
                yaml_file.write(f"{yaml_key}: {yaml_value}\n")
            elif isinstance(yaml_value, str) and yaml_value.startswith('/access/'):
                # to fit the old yaml way where paths are written without quotes.
                yaml_file.write(f"{yaml_key}: {yaml_value}\n")
            else:  # all other cases treated as strings which are single-quoted
                yaml_file.write(f"{yaml_key}: '{yaml_value}'\n")

def main():
    args = parse_args()
    # meta is a list of per-sample dicts — one dict per CSV row
    # e.g. single: [{"id": "HD829", "type": "T", ...}]
    # e.g. paired: [{"id": "HD829", "type": "T", ...}, {"id": "26MD", "type": "N", ...}]

    with open(args.meta) as meta_file:
        all_samples_meta = json.load(meta_file)

    # output_files maps Coyote YAML keys to their staged filenames
    # e.g. {"cnv": "HD829.cnvs.merged.json", "cov": "HD829.cov.json", ...}
    with open(args.json_info) as info_file:
        output_files = json.load(info_file)

    # identify tumor and normal sample from the list in all_samples_meta
    tumor_sample    = next(s for s in all_samples_meta if s['type'] in ('T', 'tumor'))
    normal_sample   = next((s for s in all_samples_meta if s['type'] in ('N', 'normal')), None)

    is_paired      = normal_sample is not None

    # for paired samples coyote expects the group name suffixed with 'p'
    coyote_group_name = args.group + 'p' if is_paired else args.group

    # reads and purity are converted to int and float types, respectively, for YAML output
    tumor_reads  = int(tumor_sample['reads'])  if tumor_sample['reads']  else None
    normal_reads = int(normal_sample['reads']) if is_paired and normal_sample['reads'] else None

    tumor_purity   = float(tumor_sample['purity']) if tumor_sample['purity'] else None
    normal_purity  = float(normal_sample['purity']) if is_paired and normal_sample['purity'] else None

    # fixed fields — always present regardless of profile
    coyote_doc = {
        'subpanel':                tumor_sample['diagnosis'],
        'name':                    coyote_group_name,
        'clarity_case_id':         tumor_sample['clarity_sample_id'],
        'clarity_control_id':      normal_sample['clarity_sample_id']  if is_paired else None,
        'clarity_case_pool_id':    tumor_sample['clarity_pool_id'],
        'clarity_control_pool_id': normal_sample['clarity_pool_id']    if is_paired else None,
        'genome_build':            38,
        'vcf_files':               build_access_path(args.subdir, 'vcf', args.vcf),
        'sample_no':               len(all_samples_meta),
        'case_id':                 tumor_sample['id'],
        'control_id':              normal_sample['id']                 if is_paired else None,
        'profile':                 args.environment,
        'assay':                   args.assay,
        'sequencing_scope':        'panel',
        'omics_layer':             'DNA',
        'sequencing_technology':   'Illumina',
        'pipeline':                args.pipeline_name,
        'pipeline_version':        args.pipeline_version,
        'case_ffpe':               tumor_sample['ffpe'],
        'case_sequencing_run':     tumor_sample['sequencing_run'],
        'case_reads':              tumor_reads,
        'case_purity':             tumor_purity,
        'control_ffpe':            normal_sample['ffpe']                if is_paired else None,
        'control_sequencing_run':  normal_sample['sequencing_run']      if is_paired else None,
        'control_reads':           normal_reads                         if is_paired else None,
        'control_purity':          normal_purity,
        'paired':                  is_paired,
    }

    # cnvprofile: gatk (modeled.png) preferred over cnvkit (cnvkit_overview.png)
    # only one cnvprofile field appears in the final YAML
    if 'cnvprofile_gatk' in output_files:
        coyote_doc['cnvprofile'] = build_access_path(args.subdir, 'plots', output_files['cnvprofile_gatk'])
    elif 'cnvprofile_cnvkit' in output_files:
        coyote_doc['cnvprofile'] = build_access_path(args.subdir, 'plots', output_files['cnvprofile_cnvkit'])

    # optional file fields — only added if present in output_files (json_INFO from OUTPUT_FILES process)
    optional_file_subfolders = {
        'cnv':        'cnv',
        'transloc':   'fusions',
        'biomarkers': 'biomarkers',
        'cov':        'QC',
        'lowcov':     'QC',
    }

    for coyote_key, subfolder in optional_file_subfolders.items():
        if coyote_key in output_files:
            coyote_doc[coyote_key] = build_access_path(args.subdir, subfolder, output_files[coyote_key])

    # purity is also written as a standalone field at the end
    # this mirrors the existing line in the COYOTE_YAML process
    if tumor_purity is not None:
        coyote_doc['purity'] = tumor_purity

    write_yaml(coyote_doc, args.out)

if __name__ == '__main__':
    main()