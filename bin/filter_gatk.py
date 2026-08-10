#!/usr/bin/env python3
"""Convert GATK copy-number records to SV-style VCF records."""

from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import add_info, leading_float, open_output, write_variant


def filter_gatk(vcf_file, out_file):
    """Convert GATK copy-number records to SV-style VCF records."""
    inserted_info_headers = False
    with vcf.VCFReader(vcf_file) as reader, open_output(out_file) as out_fh:
        for line in reader.meta_header_lines:
            if line.startswith("##FORMAT=<ID=GT"):
                print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Segment genotype">', file=out_fh)
            elif line.startswith("##INFO") and not inserted_info_headers:
                print(line, file=out_fh)
                print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', file=out_fh)
                print(
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
                    file=out_fh,
                )
                print(
                    '##INFO=<ID=gatkCN,Number=1,Type=Integer,Description="estimated copy number 0-5 gatk">',
                    file=out_fh,
                )
                inserted_info_headers = True
            else:
                print(line, file=out_fh)
        print(reader.column_header_line, file=out_fh)

        for var in reader:
            gt = var["GT"][0]
            if leading_float(gt.get("GT")) == 0:
                continue
            if leading_float(gt.get("GT")) == 1:
                gt["GT"] = "1/1" if leading_float(gt.get("CN")) == 0 else "0/1"
            svtype = var["ALT"].replace("<", "").replace(">", "")
            add_info(var, "SVTYPE", svtype)
            add_info(var, "SVLEN", str(int(leading_float(var["INFO"].get("END"))) - int(var["POS"]) + 1))
            add_info(var, "gatkCN", gt.get("CN", ""))
            write_variant(var, out_fh, use_original_csq=True)


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input GATK VCF.")
    parser.add_argument("--out", required=True, help="Output VCF.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    filter_gatk(args.vcf, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
