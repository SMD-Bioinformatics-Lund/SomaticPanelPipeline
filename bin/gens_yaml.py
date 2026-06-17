#!/usr/bin/env python3

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate a Gens (v4) somatic YAML for a case from SomaticPanelPipeline outputs'
    )
    parser.add_argument('--case_id',        required=True, help='Gens case_id (group, suffixed with p if paired)')
    parser.add_argument('--genome_build',   required=True, type=int, help='Genome build, e.g. 38')
    parser.add_argument('--gens_accessdir', required=True, help='Base /access/ directory for Gens plot data (e.g. /access/<subdir>/gens)')
    parser.add_argument('--sample_ids',     required=True, nargs='+', help='Sample IDs, one per sample')
    parser.add_argument('--sample_types',   required=True, nargs='+', help='Sample types (T/N), one per sample, same order as --sample_ids')
    parser.add_argument('--sample_sexes',   required=True, nargs='+', help='Sample sexes, one per sample, same order as --sample_ids (. if unknown)')
    parser.add_argument('--baf_filenames',  required=True, nargs='+', help='BAF bed.gz filenames, one per sample, same order as --sample_ids')
    parser.add_argument('--cov_filenames',  required=True, nargs='+', help='Coverage bed.gz filenames, one per sample, same order as --sample_ids')
    parser.add_argument('--loh_bed_filename', required=False, default=None, help='LOH category bed filename (tumor sample only); provided only for solid and GMSHem')
    parser.add_argument('--out',            required=True, help='Output YAML filename')
    return parser.parse_args()


def build_access_path(gens_accessdir, filename):
    """Build the /access/ path used by Gens to locate a plot-data file."""
    return f"{gens_accessdir}/{filename}"


def check_loh_bed_has_content(loh_bed_filename):
    """
    Returns True if the loh_bed file has content and should be referenced in the YAML.

    CALL_LOH always emits a *.loh_cat.bed file: real content when params.loh is true,
    an empty (touched) file otherwise. A missing path or a 0-byte file is
    treated as unusable.
    """
    if not loh_bed_filename:
        return False
    if not os.path.exists(loh_bed_filename):
        return False
    return os.path.getsize(loh_bed_filename) > 0


def build_sample_entries(sample_ids, sample_types, sample_sexes, baf_filenames, cov_filenames,
                          gens_accessdir, loh_bed_filename, loh_annotation_usable):
    """
    Build the list of per-sample dicts for the YAML 'samples' section.

    Each entry contains baf/coverage access paths, sample_type, and sex. For the
    tumor sample only, if loh_annotation_usable is True, a sample_annotations
    entry pointing at the LOH category file is also added.

    """
    sample_entries = []

    # zip = built-in Python function that takes multiple iterables and pairs up their elements positionally
    # 
    for sample_id, sample_type, sample_sex, baf_filename, cov_filename in zip(
        sample_ids, sample_types, sample_sexes, baf_filenames, cov_filenames
    ):
        sample_entry = {
            'sample_id':   sample_id,
            'baf':         build_access_path(gens_accessdir, baf_filename),
            'coverage':    build_access_path(gens_accessdir, cov_filename),
            'sample_type': sample_type,
            'sex':         sample_sex,
        }

        is_tumor_sample = sample_type in ('T', 'tumor')
        if loh_annotation_usable and is_tumor_sample:
            sample_entry['sample_annotations'] = [
                {
                    'file': build_access_path(gens_accessdir, loh_bed_filename),
                    'name': 'LOH',
                }
            ]

        sample_entries.append(sample_entry)

    return sample_entries


def write_gens_yaml(case_id, genome_build, sample_entries, out_filepath):
    """
    Write the nested Gens somatic YAML:

    case_id: <case_id>
    genome_build: <int>
    samples:
    - sample_id: <id>
      baf: /access/...
      coverage: /access/...
      sample_type: <T/N>
      sex: <sex>
      sample_annotations:        # optional, only for samples that have one
      - file: /access/...
        name: LOH
    """
    with open(out_filepath, 'w') as yaml_file:
        yaml_file.write('---\n')
        yaml_file.write(f"case_id: '{case_id}'\n")
        yaml_file.write(f"genome_build: {genome_build}\n")
        yaml_file.write("samples:\n")

        for sample_entry in sample_entries:
            yaml_file.write(f"- sample_id: '{sample_entry['sample_id']}'\n")
            yaml_file.write(f"  baf: {sample_entry['baf']}\n")
            yaml_file.write(f"  coverage: {sample_entry['coverage']}\n")
            yaml_file.write(f"  sample_type: '{sample_entry['sample_type']}'\n")
            yaml_file.write(f"  sex: '{sample_entry['sex']}'\n")

            sample_annotations = sample_entry.get('sample_annotations')
            if sample_annotations:
                yaml_file.write("  sample_annotations:\n")
                for annotation in sample_annotations:
                    yaml_file.write(f"  - file: {annotation['file']}\n")
                    yaml_file.write(f"    name: {annotation['name']}\n")


def main():
    args = parse_args()

    sample_count = len(args.sample_ids)
    per_sample_args = (
        ('--sample_types',  args.sample_types),
        ('--sample_sexes',  args.sample_sexes),
        ('--baf_filenames', args.baf_filenames),
        ('--cov_filenames', args.cov_filenames),
    )
    for arg_name, values in per_sample_args:
        if len(values) != sample_count:
            raise ValueError(
                f"Length mismatch: --sample_ids has {sample_count} entries, "
                f"{arg_name} has {len(values)}"
            )

    loh_annotation_usable = check_loh_bed_has_content(args.loh_bed_filename)

    sample_entries = build_sample_entries(
        sample_ids=args.sample_ids,
        sample_types=args.sample_types,
        sample_sexes=args.sample_sexes,
        baf_filenames=args.baf_filenames,
        cov_filenames=args.cov_filenames,
        gens_accessdir=args.gens_accessdir,
        loh_bed_filename=args.loh_bed_filename,
        loh_annotation_usable=loh_annotation_usable,
    )

    write_gens_yaml(args.case_id, args.genome_build, sample_entries, args.out)


if __name__ == '__main__':
    main()