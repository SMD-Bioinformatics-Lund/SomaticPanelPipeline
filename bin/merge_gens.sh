#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  merge_gens.sh --sample-id ID --case-id ID --sample-type TYPE --gens-accessdir DIR --prefix PREFIX [--split]
USAGE
}

sample_id=
case_id=
sample_type=
gens_accessdir=
prefix=
split=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-id)
            sample_id="$2"
            shift 2
            ;;
        --case-id)
            case_id="$2"
            shift 2
            ;;
        --sample-type)
            sample_type="$2"
            shift 2
            ;;
        --gens-accessdir)
            gens_accessdir="$2"
            shift 2
            ;;
        --prefix)
            prefix="$2"
            shift 2
            ;;
        --split)
            split=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 2
            ;;
    esac
done

if [[ -z "$sample_id" || -z "$case_id" || -z "$sample_type" || -z "$gens_accessdir" || -z "$prefix" ]]; then
    usage
    exit 2
fi

sort_bed() {
    local infile="$1"
    local outfile="$2"

    bedtools sort -i "$infile" > "$outfile"
}

write_base_bed() {
    local kind="$1"
    local base_bed="${prefix}.base.${kind}.bed"
    local bed_file
    local found=false

    : > "$base_bed"

    for bed_file in *."${kind}".bed.gz; do
        [[ -e "$bed_file" ]] || continue
        found=true
        zcat -f "$bed_file" | awk '/^o_/ { sub(/^o_/, ""); print }' >> "$base_bed"
    done

    if [[ "$found" == false ]]; then
        echo "No *.${kind}.bed.gz files found for split GENS merge" >&2
        exit 1
    fi
}

expand_alleles() {
    local infile="$1"
    local outfile="$2"
    local allele

    : > "$outfile"

    for allele in o a b c d; do
        awk -v allele="$allele" '{ print allele "_" $0 }' "$infile" >> "$outfile"
    done
}

merge_split_beds() {
    local kind
    local base_bed
    local sorted_base_bed
    local merged_bed
    local sorted_merged_bed

    for kind in cov baf; do
        base_bed="${prefix}.base.${kind}.bed"
        sorted_base_bed="${prefix}.base.${kind}.bed.sort"
        merged_bed="${prefix}.merged.${kind}.bed"
        sorted_merged_bed="${prefix}.merged.sorted.${kind}.bed"

        write_base_bed "$kind"
        sort_bed "$base_bed" "$sorted_base_bed"
        expand_alleles "$sorted_base_bed" "$merged_bed"
        sort_bed "$merged_bed" "$sorted_merged_bed"
        bgzip -f "$sorted_merged_bed"
        tabix -f "${sorted_merged_bed}.gz"
    done
}

index_full_beds() {
    local kind
    local bed_file

    for kind in baf cov; do
        bed_file="${prefix}.full.${kind}.bed.gz"
        if [[ ! -f "$bed_file" ]]; then
            echo "Expected GENS BED file not found: $bed_file" >&2
            exit 1
        fi
        tabix -f "$bed_file"
    done
}

write_gens_commands() {
    local filename_mod="$1"
    local baf="${gens_accessdir}/${prefix}.${filename_mod}.baf.bed.gz"
    local cov="${gens_accessdir}/${prefix}.${filename_mod}.cov.bed.gz"

    printf 'gens load sample --sample-id %s --case-id %s --genome-build 38 --baf %s --coverage %s\n' \
        "$sample_id" "$case_id" "$baf" "$cov" > "${prefix}.gens"

    printf 'gens load sample --sample-id %s --case-id %s --genome-build 38 --sample-type %s --baf %s --coverage %s\n' \
        "$sample_id" "$case_id" "$sample_type" "$baf" "$cov" > "${prefix}.gens_v4_somatic"
}

if [[ "$split" == true ]]; then
    merge_split_beds
    write_gens_commands "merged.sorted"
else
    index_full_beds
    write_gens_commands "full"
fi
