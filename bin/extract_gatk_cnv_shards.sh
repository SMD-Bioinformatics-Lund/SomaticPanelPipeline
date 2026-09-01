#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  extract_gatk_cnv_shards.sh --ploidy PLOIDY_TAR --calls-dir DIR CALL_SHARD_TAR [CALL_SHARD_TAR ...]
USAGE
}

ploidy=
calls_dir=

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ploidy)
            ploidy="$2"
            shift 2
            ;;
        --calls-dir)
            calls_dir="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown argument: $1" >&2
            usage
            exit 2
            ;;
        *)
            break
            ;;
    esac
done

if [[ -z "$ploidy" || -z "$calls_dir" || $# -eq 0 ]]; then
    usage
    exit 2
fi

tar -xf "$ploidy"
mkdir -p "$calls_dir"

for shard_tar in "$@"; do
    tar -xf "$shard_tar"
done

find . -mindepth 1 -maxdepth 1 -type d -name '*_[0-9]*' -exec mv -t "$calls_dir" {} +
