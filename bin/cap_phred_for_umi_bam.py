#!/usr/bin/env python3

import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(
        description="Cap BAM base qualities at a specified Phred score."
    )
    parser.add_argument("input_bam")
    parser.add_argument("output_bam")
    parser.add_argument(
        "--max-quality",
        type=int,
        default=70,
        help="Maximum permitted base quality (default: 70)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="BAM compression threads (default: 4)",
    )
    args = parser.parse_args()

    changed_reads = 0
    changed_bases = 0

    with pysam.AlignmentFile(
        args.input_bam, "rb", threads=args.threads
    ) as infile, pysam.AlignmentFile(
        args.output_bam,
        "wb",
        template=infile,
        threads=args.threads,
    ) as outfile:

        for read in infile.fetch(until_eof=True):
            qualities = read.query_qualities

            # QUAL="*" is represented as None; leave it unchanged.
            if qualities is not None:
                capped = [
                    min(q, args.max_quality)
                    for q in qualities
                ]

                number_changed = sum(
                    q > args.max_quality for q in qualities
                )

                if number_changed:
                    read.query_qualities = capped
                    changed_reads += 1
                    changed_bases += number_changed

            outfile.write(read)

    print(f"Reads modified: {changed_reads:,}")
    print(f"Base qualities capped: {changed_bases:,}")


if __name__ == "__main__":
    main()
