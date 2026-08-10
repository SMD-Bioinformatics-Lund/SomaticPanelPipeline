#!/usr/bin/env python3
"""Parse sambamba depth output and report contiguous low-depth regions."""

from argparse import ArgumentParser
from pathlib import Path

from vcf_pipeline_utils import leading_float, open_output


def panel_depth(depth_file, out_file, cutoff=500):
    """Report contiguous low-depth panel regions from sambamba depth base output."""
    if not Path(depth_file).is_file():
        raise FileNotFoundError(f"No file {depth_file}")

    with open(depth_file) as in_fh, open_output(out_file) as out_fh:
        start_pos = None
        start_chr = None
        last_low_pos = None
        last_low_chr = None
        low_cov_sum = 0.0

        for line in in_fh:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                continue

            chrom, pos_text, depth_text = fields[:3]
            try:
                pos = int(pos_text)
            except ValueError:
                continue

            depth = leading_float(depth_text)

            if depth < cutoff:
                if last_low_chr == chrom and last_low_pos == pos - 1:
                    low_cov_sum += depth
                else:
                    if start_pos is not None and start_chr is not None:
                        out_fh.write(
                            f"{start_chr}\t{start_pos}\t{last_low_pos}\t"
                            f"{low_cov_sum / (last_low_pos - start_pos + 1)}\n"
                        )
                    start_chr = chrom
                    start_pos = pos
                    low_cov_sum = depth

                last_low_chr = chrom
                last_low_pos = pos

            elif start_chr is not None and start_pos is not None:
                out_fh.write(
                    f"{start_chr}\t{start_pos}\t{last_low_pos}\t"
                    f"{low_cov_sum / (last_low_pos - start_pos + 1)}\n"
                )
                start_chr = None
                start_pos = None
                last_low_chr = None
                last_low_pos = None
                low_cov_sum = 0.0

        if start_chr is not None and start_pos is not None:
            out_fh.write(
                f"{start_chr}\t{start_pos}\t{last_low_pos}\t"
                f"{low_cov_sum / (last_low_pos - start_pos + 1)}\n"
            )


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True, help="Input TSV from sambamba depth base.")
    parser.add_argument("--cutoff", type=float, default=500, help="Low-depth cutoff.")
    parser.add_argument("--out", required=True, help="Output BED-like low-depth regions.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    panel_depth(args.depth, args.out, args.cutoff)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
