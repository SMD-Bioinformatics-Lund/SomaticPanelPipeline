#!/usr/bin/env python3
"""Create Coyote CN segment BED files from merged SV calls."""

from argparse import ArgumentParser

import vcf2 as vcf
from vcf_pipeline_utils import leading_float


def _read_panel(path, panel_genes):
    """Read one coyote segmentator panel file into a shared dictionary."""
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                panel_genes[cols[0]] = {"TYPE": cols[1], "QUESTION": cols[2]}


def _panel_gene_hits(genes, panel_genes):
    """Return matching panel gene annotations."""
    interesting = []
    for gene in genes.split(","):
        if gene in panel_genes:
            interesting.append(f"{gene}:{panel_genes[gene]['TYPE']}:{panel_genes[gene]['QUESTION']}")
    return len(interesting), ",".join(interesting)


def _header_callers(header_text):
    """Infer caller names from a merged SV VCF header."""
    callers = []
    for caller in ("gatk", "manta", "cnvkit"):
        if caller in header_text:
            callers.append(caller)
    return callers


def _manta_pr_sr(var, sample_id):
    """Return Manta PR and SR strings for one sample."""
    for gt in var["GT"]:
        if gt["_sample_id"] == sample_id:
            return gt.get("PR", "").replace(",", "/"), gt.get("SR", "").replace(",", "/")
    return 0, 0


def write_coyote_segment_intervals(vcf_file, sample_id, out_file):
    """Write unannotated Coyote CN segment intervals from a merged SV VCF."""
    with vcf.VCFReader(vcf_file) as reader, open(out_file, "w") as out_fh:
        callers_from_header = _header_callers(reader.header_str)
        for var in reader:
            svtype = var["INFO"].get("SVTYPE")
            if svtype == "BND":
                continue
            start = int(var["POS"])
            if var["INFO"].get("END"):
                end = int(leading_float(var["INFO"]["END"]))
            else:
                end = start + abs(int(leading_float(var["INFO"].get("SVLEN"))))
            probes = var["INFO"].get("PROBES", 0)
            if var["INFO"].get("FOLD_CHANGE_LOG"):
                fold = var["INFO"]["FOLD_CHANGE_LOG"]
            elif var["INFO"].get("gatkCN"):
                fold = leading_float(var["INFO"]["gatkCN"]) / 2
            elif svtype == "DEL":
                fold = "DEL"
            else:
                fold = "AMP"
            callers = (
                "&".join(callers_from_header)
                if var["INFO"].get("set") == "Intersection"
                else var["INFO"].get("set", "")
            )
            sign = "-" if svtype == "DEL" else "+"
            pr, sr = _manta_pr_sr(var, sample_id) if "manta" in callers else (0, 0)
            print(
                f"{var['CHROM']}\t{start}\t{end}\t{probes}\t{fold}\t{sign}\t{callers}\t{pr}:{sr}",
                file=out_fh,
            )
    return out_file


def coyote_segmentator(
    intersect_bed,
    sample_id,
    panel_files,
    normal=False,
    out_prefix=None,
    panel_out=None,
    raw_out=None,
):
    """Create coyote CN segment BED files from precomputed bedtools intersect output."""
    panel_genes = {}
    for panel in panel_files.split(","):
        _read_panel(panel, panel_genes)

    out_prefix = out_prefix or sample_id
    filtered_path = panel_out or f"{out_prefix}.cn-segments.panel.bed"
    unfiltered_path = raw_out or f"{out_prefix}.cn-segments.bed"
    groups = []
    previous_region = None
    genes = []
    previous_bed = None
    with open(intersect_bed) as in_fh:
        for line in in_fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            region = f"{fields[0]}\t{fields[1]}\t{fields[2]}\t"
            bed_str = fields[:8]
            if previous_region and region != previous_region:
                groups.append((previous_bed, genes))
                genes = []
            genes.append(fields[-1])
            previous_region = region
            previous_bed = bed_str
    if previous_bed:
        groups.append((previous_bed, genes))

    with open(filtered_path, "w") as filtered_fh, open(unfiltered_path, "w") as unfiltered_fh:
        for bed_fields, gene_list in groups:
            gene_conc = ",".join(gene_list)
            match_count, interesting = _panel_gene_hits(gene_conc, panel_genes)
            if match_count > 0:
                filtered_fh.write("\t".join(bed_fields[:6]))
                filtered_fh.write(f"\t{gene_conc}\t{interesting}\t")
                filtered_fh.write("\t".join(bed_fields[6:8]))
                if normal:
                    filtered_fh.write("\tNORMAL")
                filtered_fh.write("\n")
            unfiltered_fh.write("\t".join(bed_fields[:6]))
            unfiltered_fh.write(f"\t{gene_conc}\t")
            unfiltered_fh.write("\t".join(bed_fields[6:8]) + "\n")
    return filtered_path, unfiltered_path


def build_parser():
    parser = ArgumentParser(description="Create coyote CN segment BED files from a merged SV VCF.")
    parser.add_argument("--vcf", help="Merged SV VCF.")
    parser.add_argument("--segments-out", help="Write unannotated segment intervals and exit.")
    parser.add_argument("--intersect-bed", help="Input BED produced by bedtools intersect -loj.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--panel", help="Comma-delimited panel definition files.")
    parser.add_argument("--normal", action="store_true", help="Mark filtered segments as NORMAL.")
    parser.add_argument("--out-prefix", help="Output prefix. Defaults to sample ID.")
    parser.add_argument("--panel-out", help="Output BED with panel-matched segments.")
    parser.add_argument("--raw-out", help="Output BED with all segments.")
    return parser

def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.segments_out:
        if not args.vcf:
            raise SystemExit("--vcf is required with --segments-out")
        write_coyote_segment_intervals(args.vcf, args.sample_id, args.segments_out)
        return 0
    if not args.intersect_bed:
        raise SystemExit("--intersect-bed is required unless --segments-out is used")
    if not args.panel:
        raise SystemExit("--panel is required when formatting intersected segments")
    coyote_segmentator(
        args.intersect_bed,
        args.sample_id,
        args.panel,
        args.normal,
        args.out_prefix,
        args.panel_out,
        args.raw_out,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
