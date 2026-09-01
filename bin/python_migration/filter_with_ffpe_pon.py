#!/usr/bin/env python3
"""Annotate and filter variants using FFPE panel-of-normals files."""

from argparse import ArgumentParser

from filter_with_pon import parse_pons, read_pon  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf2  # noqa: E402
from vcf_pipeline_utils import add_info, leading_float, open_output, write_filter_header, write_variant  # noqa: E402


def compact_count(value):
    """Format a count value like Perl scalar arithmetic output."""
    return f"{value:g}"


def pon_header_lines(callers):
    """Build FFPE PON header lines.

    Args:
        callers: Iterable of caller names.

    Returns:
        List of VCF header lines.
    """
    lines = []
    for caller in callers:
        lines.extend(
            [
                f'##FILTER=<ID=WARN_FFPE_PON_{caller},Description="Warning, exists in PON for {caller}">',
                f'##FILTER=<ID=FAIL_FFPE_PON_{caller},Description="Failed, exists in PON for {caller}">',
                (
                    f'##INFO=<ID=PON_FFPE_NUM_{caller},Number=.,Type=String,'
                    f'Description="Number of PON samples variants exists in, for {caller}">'
                ),
                f'##INFO=<ID=PON_FFPE_VAFS_{caller},Number=.,Type=String,Description="VAFs for PON variants, for {caller}">',
            ]
        )
    return lines


def filter_vcf(vcf_file, tumor_id, pons_arg, out_file):
    """Filter a VCF using FFPE panel-of-normals annotations."""
    pon_files = parse_pons(pons_arg)
    pons = {caller: read_pon(path) for caller, path in pon_files.items()}

    with vcf2.VCFReader(file=vcf_file) as reader, open_output(out_file) as out_fh:
        write_filter_header(reader, out_fh, pon_header_lines(pons.keys()))
        for var in reader:
            tumor_vaf = None
            for gt in var["GT"]:
                if tumor_id and gt["_sample_id"] == tumor_id:
                    tumor_vaf = leading_float(gt.get("VAF"))
            if tumor_vaf is None:
                raise ValueError("Tumor sample not found")

            varid = f"{var['CHROM']}_{var['POS']}_{var['REF']}_{var['ALT']}"
            filters = var["FILTER"].split(";")
            remove_pass = False
            for caller in sorted(pons):
                pon_item = pons[caller].get(varid)
                if not pon_item:
                    continue
                frac = pon_item["num_non_germ"] / pon_item["total"]
                if frac > 0.05:
                    if tumor_vaf < pon_item["mean_vaf"] + pon_item["stdev_vaf"] * 2:
                        filters.append(f"FAIL_FFPE_PON_{caller}")
                        remove_pass = True
                    else:
                        filters.append(f"WARN_FFPE_PON_{caller}")
                add_info(
                    var,
                    f"PON_FFPE_NUM_{caller}",
                    f"{compact_count(pon_item['num_non_germ'])}/{compact_count(pon_item['total'])}",
                )
                add_info(var, f"PON_FFPE_VAFS_{caller}", pon_item["all"])

            if remove_pass:
                filters = [item for item in filters if item != "PASS"]
            var["FILTER"] = ";".join(filters)
            write_variant(var, out_fh, use_original_csq=True, trailing_sample_tab=False)


def build_parser():
    """Build the command-line parser."""
    parser = ArgumentParser(description="Filter variants with FFPE panel-of-normals files.")
    parser.add_argument("--vcf", required=True, help="Input VCF file.")
    parser.add_argument("--tumor-id", required=True, help="Tumor sample ID.")
    parser.add_argument("--pons", required=True, help="Comma-delimited caller=PON table entries.")
    parser.add_argument("--out", required=True, help="Output VCF file.")
    return parser


def main(argv=None):
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    filter_vcf(args.vcf, args.tumor_id, args.pons, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
