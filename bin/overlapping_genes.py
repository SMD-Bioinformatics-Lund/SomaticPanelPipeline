#!/usr/bin/env python3
"""Parse bedtools-intersected low-coverage BED regions and annotate with gene names."""

from __future__ import annotations

from argparse import ArgumentParser

from vcf_pipeline_utils import open_output


def overlapping_genes(intersect_bed, out_file):
    """Parse precomputed bedtools intersect output and write lowcov regions with genes.

    Expected input is already produced by something like:

        bedtools intersect -a lowcov.bed -b genes.bed -loj > lowcov.intersect.bed

    Expected lowcov.bed columns:

        chrom   start   end   mean_depth

    Expected bedtools -loj output contains the original lowcov columns first,
    followed by the overlapping gene BED columns.

    Output:

        chrom   start   end   mean_depth   gene
    """
    with open(intersect_bed) as in_fh, open_output(out_file) as out_fh:
        for line in in_fh:
            line = line.rstrip("\n")

            if not line:
                continue

            if (
                line.startswith("@")
                or line.startswith("CONTIG")
                or line.startswith("REF")
            ):
                continue

            fields = line.split("\t")

            # Need at least the original lowcov BED columns.
            if len(fields) < 4:
                continue

            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            mean_depth = fields[3]

            # With bedtools -loj, no overlap usually gives "." / -1 fields
            # at the end. The old script used fields[-1] as the gene name,
            # so we keep the same behavior.
            gene = fields[-1] if len(fields) > 4 else "."

            out_fh.write(f"{chrom}\t{start}\t{end}\t{mean_depth}\t{gene}\n")


def build_parser():
    parser = ArgumentParser(
        description="Parse precomputed bedtools-intersected lowcov BED and annotate with overlapping genes."
    )
    parser.add_argument(
        "--intersect-bed",
        required=True,
        help="Input BED produced by bedtools intersect -loj.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output lowcov BED with gene annotation.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    overlapping_genes(args.intersect_bed, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
