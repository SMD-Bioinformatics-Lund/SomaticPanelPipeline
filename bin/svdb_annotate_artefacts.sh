#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  svdb_annotate_artefacts.sh --query-vcf VCF --dbs DBS_TSV --out OUT -- [SVDB_ARGS...]

DBS_TSV must contain two tab-separated columns:
  database_name    database_vcf
USAGE
}

query_vcf=
dbs_tsv=
out_vcf=

while [[ $# -gt 0 ]]; do
    case "$1" in
        --query-vcf)
            query_vcf="$2"
            shift 2
            ;;
        --dbs)
            dbs_tsv="$2"
            shift 2
            ;;
        --out)
            out_vcf="$2"
            shift 2
            ;;
        --)
            shift
            break
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

if [[ -z "$query_vcf" || -z "$dbs_tsv" || -z "$out_vcf" ]]; then
    usage
    exit 2
fi

input_vcf="$query_vcf"

while IFS=$'\t' read -r db_name db_vcf; do
    [[ -n "$db_name" ]] || continue

    tmp_vcf="tmp_${db_name}.vcf"
    svdb --query "$@" \
        --out_frq "AFRQ_${db_name}" \
        --out_occ "ACOUNT_${db_name}" \
        --db "$db_vcf" \
        --query_vcf "$input_vcf" \
        > "$tmp_vcf"

    input_vcf="$tmp_vcf"
done < "$dbs_tsv"

mv "$input_vcf" "$out_vcf"
