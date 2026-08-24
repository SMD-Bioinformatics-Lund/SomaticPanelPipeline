#!/usr/bin/env python3
"""
Description: this script generates a YAML document for Coyote3 data import, combining meta info, optional file locations and other
parameters. This script was created to replace the currently used processes "COYOTE" and "COYOTE_YAML".
"""

import argparse
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a COYOTE YAML document from SomaticPanelPipeline outputs and all_samples_metadata"
    )
    parser.add_argument("--group", required=True, help="sample group name")
    parser.add_argument(
        "--meta", required=True, help="path to meta JSON file - one dict per sample"
    )
    parser.add_argument("--vcf", required=True, help="VCF filename")
    parser.add_argument(
        "--json_info",
        required=True,
        help="path to JSON file containing output filenames",
    )
    parser.add_argument(
        "--subdir",
        required=True,
        help="pipeline subdirectory for building /access/ paths",
    )
    parser.add_argument("--assay", required=True, help="coyote group/assay name")
    parser.add_argument(
        "--environment",
        required=True,
        help="production, development, validation, or testing",
    )
    parser.add_argument(
        "--pipeline_name", required=True, help="pipeline name from workflow manifest"
    )
    parser.add_argument(
        "--pipeline_version",
        required=True,
        help="pipeline version from workflow manifest",
    )
    parser.add_argument("--out", required=True, help="output YAML filename")
    return parser.parse_args()


def build_access_path(subdir, subfolder, filename):
    """Build the /access/ path used by Coyote to locate a file."""
    return f"/access/{subdir}/{subfolder}/{filename}"


def sample_type(sample):
    return str(sample.get("type", "")).lower()


def optional_int(value):
    if value in (None, "", False):
        return None
    if isinstance(value, str) and value.lower() == "false":
        return None
    return int(value)


def optional_float(value):
    if value in (None, "", False):
        return None
    if isinstance(value, str) and value.lower() == "false":
        return None
    return float(value)


def normalize_sex(value):
    sex = str(value or "").strip().lower()
    if sex in ("male", "m"):
        return "male"
    if sex in ("female", "f"):
        return "female"
    return "unknown"


def write_yaml(coyote_doc, out_filepath):
    """
    Write a flat YAML file.
    Handles None, bool, int/float and string values with correct YAML formatting.
    """
    with open(out_filepath, "w", encoding="utf-8") as yaml_file:
        yaml_file.write("---\n")  # Start of YAML document
        for yaml_key, yaml_value in coyote_doc.items():
            if yaml_value is None:
                yaml_file.write(f"{yaml_key}: null\n")
            elif isinstance(yaml_value, bool):
                yaml_file.write(f"{yaml_key}: {'true' if yaml_value else 'false'}\n")
            elif isinstance(yaml_value, (int, float)):
                yaml_file.write(f"{yaml_key}: {yaml_value}\n")
            elif isinstance(yaml_value, str) and yaml_value.startswith("/access/"):
                yaml_file.write(f"{yaml_key}: {yaml_value}\n")
            else:
                yaml_file.write(f"{yaml_key}: '{yaml_value}'\n")


def main():
    args = parse_args()

    with open(args.meta, encoding="utf-8") as meta_file:
        all_samples_meta = json.load(meta_file)

    with open(args.json_info, encoding="utf-8") as info_file:
        output_files = json.load(info_file)

    tumor_sample = next(
        (s for s in all_samples_meta if sample_type(s) in ("t", "tumor")),
        all_samples_meta[0],
    )
    normal_sample = next(
        (s for s in all_samples_meta if sample_type(s) in ("n", "normal")), None
    )

    is_paired = normal_sample is not None

    coyote_group_name = args.group + "p" if is_paired else args.group

    tumor_reads = optional_int(tumor_sample.get("reads"))
    normal_reads = optional_int(normal_sample.get("reads")) if is_paired else None

    tumor_purity = optional_float(tumor_sample.get("purity"))
    normal_purity = optional_float(normal_sample.get("purity")) if is_paired else None

    coyote_doc = {
        "subpanel": tumor_sample["diagnosis"],
        "name": coyote_group_name,
        "clarity_case_id": tumor_sample["clarity_sample_id"],
        "clarity_control_id": normal_sample["clarity_sample_id"] if is_paired else None,
        "clarity_case_pool_id": tumor_sample["clarity_pool_id"],
        "clarity_control_pool_id": (
            normal_sample["clarity_pool_id"] if is_paired else None
        ),
        "genome_build": 38,
        "vcf_files": build_access_path(args.subdir, "vcf", args.vcf),
        "sample_no": len(all_samples_meta),
        "case_id": tumor_sample["id"],
        "control_id": normal_sample["id"] if is_paired else None,
        "profile": args.environment,
        "assay": args.assay,
        "sequencing_scope": "panel",
        "omics_layer": "DNA",
        "sequencing_technology": "Illumina",
        "sex": normalize_sex(tumor_sample.get("sex")),
        "readmode": "PE",
        "pipeline": args.pipeline_name,
        "pipeline_version": args.pipeline_version,
        "case_ffpe": tumor_sample["ffpe"],
        "case_sequencing_run": tumor_sample["sequencing_run"],
        "case_reads": tumor_reads,
        "case_purity": tumor_purity,
        "control_ffpe": normal_sample["ffpe"] if is_paired else None,
        "control_sequencing_run": (
            normal_sample["sequencing_run"] if is_paired else None
        ),
        "control_reads": normal_reads if is_paired else None,
        "control_purity": normal_purity,
        "paired": is_paired,
    }

    if "cnvprofile_gatk" in output_files:
        coyote_doc["cnvprofile"] = build_access_path(
            args.subdir, "plots", output_files["cnvprofile_gatk"]
        )
    elif "cnvprofile_cnvkit" in output_files:
        coyote_doc["cnvprofile"] = build_access_path(
            args.subdir, "plots", output_files["cnvprofile_cnvkit"]
        )

    optional_file_subfolders = {
        "cnv": "cnv",
        "transloc": "fusions",
        "biomarkers": "biomarkers",
        "cov": "QC",
        "lowcov": "QC",
    }

    for coyote_key, subfolder in optional_file_subfolders.items():
        if coyote_key in output_files:
            coyote_doc[coyote_key] = build_access_path(
                args.subdir, subfolder, output_files[coyote_key]
            )

    if tumor_purity is not None:
        coyote_doc["purity"] = tumor_purity

    write_yaml(coyote_doc, args.out)


if __name__ == "__main__":
    main()
