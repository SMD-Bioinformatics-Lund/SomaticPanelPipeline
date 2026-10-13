#!/usr/bin/env python3
"""Transform GENS BED rows and write import commands without running external tools."""

import argparse
from pathlib import Path


def prepare_bed(args):
    with args.out.open('w', newline='') as output:
        for path in sorted(args.beds, key=lambda path: path.name):
            with path.open(newline='') as source:
                for line in source:
                    if line.startswith('o_'):
                        output.write(line[2:].removesuffix('\n') + '\n')


def expand_alleles(args):
    with args.out.open('w', newline='') as output:
        for allele in ('o', 'a', 'b', 'c', 'd'):
            with args.bed.open(newline='') as source:
                for line in source:
                    row = line.removesuffix('\n')
                    output.write(f"{allele}_{row}\n")


def write_commands(args):
    command = (
        f"gens load sample --sample-id {args.sample_id} --case-id {args.case_id} "
        f"--genome-build {args.genome_build}"
    )
    files = f" --baf {args.access_dir}/{args.baf} --coverage {args.access_dir}/{args.coverage}\n"
    Path(f'{args.prefix}.gens').write_text(command + files)
    Path(f'{args.prefix}.gens_v4_somatic').write_text(
        command + f" --sample-type {args.sample_type}" + files
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare = commands.add_parser('prepare', help='Combine original-resolution rows without changing numeric values')
    prepare.add_argument('--beds', type=Path, nargs='+', required=True)
    prepare.add_argument('--out', type=Path, required=True)
    prepare.set_defaults(run=prepare_bed)
    expand = commands.add_parser('expand', help='Add the five legacy allele prefixes')
    expand.add_argument('--bed', type=Path, required=True)
    expand.add_argument('--out', type=Path, required=True)
    expand.set_defaults(run=expand_alleles)
    writer = commands.add_parser('commands', help='Write GENS import command files')
    for name in ('sample-id', 'case-id', 'sample-type', 'access-dir', 'baf', 'coverage', 'prefix'):
        writer.add_argument(f'--{name}', required=True)
    writer.add_argument('--genome-build', type=int, default=38)
    writer.set_defaults(run=write_commands)
    args = parser.parse_args()
    args.run(args)


if __name__ == '__main__':
    main()
