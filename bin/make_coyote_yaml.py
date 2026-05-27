#!/usr/bin/env python3

import argparse
import json
import yaml

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--group',            required=True)
    parser.add_argument('--meta',             required=True, help='path to meta JSON file')
    parser.add_argument('--vcf',              required=True)
    parser.add_argument('--json_info',        required=True)
    parser.add_argument('--subdir',           required=True)
    parser.add_argument('--assay',            required=True)
    parser.add_argument('--environment',      required=True)
    parser.add_argument('--pipeline_name',    required=True)
    parser.add_argument('--pipeline_version', required=True)
    parser.add_argument('--out',              required=True)
    return parser.parse_args()

def build_access_path(subdir, subfolder, filename):
    return f"/access/{subdir}/{subfolder}/{filename}"

def main():
    args = parse_args()

    with open(args.meta) as meta_file:
        meta = json.load(meta_file)

    with open(args.json_info) as info_file:
        output_files = json.load(info_file)

    # determine tumor and normal sample indices
    sample_types = meta['type'] if isinstance(meta['type'], list) else [meta['type']]
    tumor_idx    = next(i for i, t in enumerate(sample_types) if t in ('T', 'tumor'))
    normal_idx   = next((i for i, t in enumerate(sample_types) if t in ('N', 'normal')), None)

    is_paired      = normal_idx is not None
    coyote_name    = args.group + 'p' if is_paired else args.group
    tumor_purity   = float(meta['purity'][tumor_idx]) if meta['purity'][tumor_idx] else None
    normal_purity  = float(meta['purity'][normal_idx]) if is_paired and meta['purity'][normal_idx] else None

    # fixed fields — always present regardless of profile
    coyote_doc = {
        'subpanel':                meta['diagnosis'][tumor_idx],
        'name':                    coyote_name,
        'clarity_case_id':         meta['clarity_sample_id'][tumor_idx],
        'clarity_control_id':      meta['clarity_sample_id'][normal_idx]  if is_paired else None,
        'clarity_case_pool_id':    meta['clarity_pool_id'][tumor_idx],
        'clarity_control_pool_id': meta['clarity_pool_id'][normal_idx]    if is_paired else None,
        'genome_build':            38,
        'vcf_files':               build_access_path(args.subdir, 'vcf', args.vcf),
        'sample_no':               len(meta['id']),
        'case_id':                 meta['id'][tumor_idx],
        'control_id':              meta['id'][normal_idx]                  if is_paired else None,
        'profile':                 args.environment,
        'assay':                   args.assay,
        'sequencing_scope':        'panel',
        'omics_layer':             'DNA',
        'sequencing_technology':   'Illumina',
        'pipeline':                args.pipeline_name,
        'pipeline_version':        args.pipeline_version,
        'case_ffpe':               meta['ffpe'][tumor_idx],
        'case_sequencing_run':     meta['sequencing_run'][tumor_idx],
        'case_reads':              meta['reads'][tumor_idx],
        'case_purity':             tumor_purity,
        'control_ffpe':            meta['ffpe'][normal_idx]                if is_paired else None,
        'control_sequencing_run':  meta['sequencing_run'][normal_idx]      if is_paired else None,
        'control_reads':           meta['reads'][normal_idx]               if is_paired else None,
        'control_purity':          normal_purity,
        'paired':                  is_paired,
    }

    # cnvprofile: gatk (modeled.png) preferred over cnvkit (cnvkit_overview.png)
    if 'cnvprofile_gatk' in output_files:
        coyote_doc['cnvprofile'] = build_access_path(args.subdir, 'plots', output_files['cnvprofile_gatk'])
    elif 'cnvprofile_cnvkit' in output_files:
        coyote_doc['cnvprofile'] = build_access_path(args.subdir, 'plots', output_files['cnvprofile_cnvkit'])

    # optional file fields — only added if present in output_files
    optional_subfolders = {
        'cnv':        'cnv',
        'transloc':   'fusions',
        'biomarkers': 'biomarkers',
        'cov':        'QC',
        'lowcov':     'QC',
    }

    for coyote_key, subfolder in optional_subfolders.items():
        if coyote_key in output_files:
            coyote_doc[coyote_key] = build_access_path(args.subdir, subfolder, output_files[coyote_key])

    # purity as standalone field — mirrors existing behaviour in COYOTE_YAML process
    if tumor_purity is not None:
        coyote_doc['purity'] = tumor_purity

    with open(args.out, 'w') as output_yaml:
        yaml.dump(coyote_doc, output_yaml,
                  default_flow_style=False,
                  allow_unicode=True,
                  sort_keys=False,
                  explicit_start=True)

if __name__ == '__main__':
    main()