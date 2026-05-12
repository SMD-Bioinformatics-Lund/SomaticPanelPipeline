"""Structural-variant filtering and conversion utilities."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from argparse import ArgumentParser
from pathlib import Path

import vcf2 as vcf
from vcf_pipeline_utils import add_info, leading_float, open_output, write_variant


def _parse_info_column(info_text):
    """Parse a raw INFO column into a dictionary."""
    info = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "defined"
    return info


def _natural_key(value):
    """Return a sort key close to GNU sort ``-V`` for chromosome strings."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


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
                    '##INFO=<ID=gatkCN,Number=1,Type=Integer,Description="estimated copy number 0-5 gatk">', file=out_fh
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


def _manta_passes_af(var, sample_id, min_af):
    """Return whether a Manta record has enough PR/SR support."""
    for gt in var["GT"]:
        if gt["_sample_id"] != sample_id:
            continue
        pr = [0.0, 0.0]
        sr = [0.0, 0.0]
        if gt.get("PR"):
            pr = [leading_float(item) for item in gt["PR"].split(",")[:2]]
        if gt.get("SR"):
            sr = [leading_float(item) for item in gt["SR"].split(",")[:2]]
        if (pr[0] == 0 and sr[0] == 0) or ((pr[0] + pr[1]) < 20 and (sr[0] + sr[1]) < 20):
            return False
        af_pr = pr[1] / (pr[0] + pr[1]) if pr[0] + pr[1] else 0
        af_sr = sr[1] / (sr[0] + sr[1]) if sr[0] + sr[1] else 0
        if af_sr < min_af and af_pr < min_af:
            return False
    return True


def filter_manta(vcf_file, sample_id, min_af, out_prefix):
    """Split and filter Manta BND and non-BND records."""
    bnd_path = f"{out_prefix}_manta_bnd_filtered.vcf"
    filtered_path = f"{out_prefix}_manta_filtered.vcf"
    with vcf.VCFReader(vcf_file) as reader, open(bnd_path, "w") as bnd_fh, open(filtered_path, "w") as filtered_fh:
        for line in reader.meta_header_lines:
            if line.startswith("##contig") and not re.search(r"ID=[0-9XYM]{1,2},", line):
                continue
            print(line, file=bnd_fh)
            print(line, file=filtered_fh)
        print(reader.column_header_line, file=bnd_fh)
        print(reader.column_header_line, file=filtered_fh)

        for var in reader:
            if not _manta_passes_af(var, sample_id, leading_float(min_af)):
                continue
            if var["INFO"].get("SVTYPE") == "BND":
                print(var["vcf_str"], file=bnd_fh)
            else:
                print(var["vcf_str"], file=filtered_fh)
    return bnd_path, filtered_path


def _format_merged_gatk_record(agg, info, tumor=False):
    """Serialize one merged GATK/CNV segment record."""
    fields = agg[:7]
    if tumor:
        info_items = []
        for key in ("END", "SVLEN", "SVTYPE", "FOLD_CHANGE_LOG", "FOLD_CHANGE", "PROBES"):
            values = str(info[key]).split(",")
            if key in {"FOLD_CHANGE_LOG", "FOLD_CHANGE"}:
                value = sum(leading_float(item) for item in values) / len(values)
            elif key == "PROBES":
                value = sum(leading_float(item) for item in values)
            else:
                value = info[key]
            info_items.append(f"{key}={value}")
    else:
        info_items = [f"{key}={info[key]}" for key in ("END", "SVLEN", "SVTYPE", "gatkCN")]
    return "\t".join(fields + [";".join(info_items)] + agg[8:])


def _merge_gatk_type(vcf_file, svtype, tumor=False):
    """Merge adjacent GATK/CNV records for one SV type."""
    headers = []
    variants = []
    agg = []
    agg_info = {}
    merge_key = "FOLD_CHANGE" if tumor else "gatkCN"
    with vcf.open_vcf(vcf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                if not tumor and re.search(r"<ID=gatkCN,", line):
                    line = line.replace("Number=1,", "Number=.,")
                headers.append(line)
                continue
            data = line.rstrip("\n").split("\t")
            info = _parse_info_column(data[7])
            if info.get("SVTYPE") != svtype:
                continue
            if agg:
                dist = abs(int(data[1]) - int(leading_float(agg_info["END"])))
                agg_first = str(agg_info[merge_key]).split(",")[0]
                if (
                    data[0] == agg[0]
                    and dist / abs(leading_float(info["SVLEN"])) < 0.1
                    and dist / abs(leading_float(agg_info["SVLEN"])) < 0.1
                    and dist < 50000
                    and abs(leading_float(agg_first) - leading_float(info[merge_key])) < 0.2
                    and info["SVTYPE"] == agg_info["SVTYPE"]
                ):
                    agg[2] += "," + data[2]
                    agg_info["END"] = info["END"]
                    agg_info["SVLEN"] = leading_float(agg_info["SVLEN"]) + leading_float(info["SVLEN"])
                    if tumor:
                        for key in ("FOLD_CHANGE_LOG", "FOLD_CHANGE", "PROBES"):
                            agg_info[key] = f"{agg_info[key]},{info[key]}"
                    else:
                        agg_info["gatkCN"] = f"{agg_info['gatkCN']},{info['gatkCN']}"
                    continue
                variants.append(_format_merged_gatk_record(agg, agg_info, tumor=tumor))
            agg = data
            agg_info = info
    if agg:
        variants.append(_format_merged_gatk_record(agg, agg_info, tumor=tumor))
    return headers, variants


def merge_gatk(vcf_file, out_file, tumor=False):
    """Merge adjacent GATK segment records."""
    headers, deletions = _merge_gatk_type(vcf_file, "DEL", tumor=tumor)
    _headers, duplications = _merge_gatk_type(vcf_file, "DUP", tumor=tumor)
    records = sorted(
        deletions + duplications, key=lambda row: (_natural_key(row.split("\t", 1)[0]), int(row.split("\t", 2)[1]))
    )
    with open_output(out_file) as out_fh:
        out_fh.writelines(headers)
        for row in records:
            out_fh.write(row if row.endswith("\n") else row + "\n")


def merge_melt(vcf_header, sample_id, out_file):
    """Merge MELT ALU, LINE1, and SVA VCF files into one Scout-style VCF."""
    qcval = {
        "lc": "Low Complex Region",
        "s25": "Greater than 25.0% of samples do not have data",
        "hDP": '"More than the expected number of discordant pairs at this site are also split',
        "PASS": "PASS",
        "rSD": "Ratio of LP to RP is greater than 2.0 standard deviations",
    }
    assess = {
        "0": "No overlapping reads at site",
        "1": "Imprecise breakpoint due to greater than expected distance between evidence",
        "2": "discordant pair evidence only -- No split read information",
        "3": "left side TSD evidence only",
        "4": "right side TSD evidence only",
        "5": "TSD decided with split reads -> highest possible quality",
    }
    mod_files = []
    for path in ("ALU.final_comp.vcf", "LINE1.final_comp.vcf", "SVA.final_comp.vcf"):
        if not Path(path).is_file() or Path(path).stat().st_size < 1:
            continue
        out_path = f"{path}mod"
        mod_files.append(out_path)
        with vcf.VCFReader(path) as reader, open(out_path, "w") as out_fh:
            inserted = False
            for line in reader.meta_header_lines:
                line = line.replace("'\\|'", "pipe-sign")
                if line.startswith("##INFO") and not inserted:
                    print(
                        '##INFO=<ID=END,Number=.,Type=Integer,Description="END position set to start position for insertions">',
                        file=out_fh,
                    )
                    print(
                        '##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description="Custom annotations for scout">',
                        file=out_fh,
                    )
                    inserted = True
                print(line, file=out_fh)
            print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}", file=out_fh)
            for var in reader:
                if "ac0" in var["FILTER"]:
                    continue
                alt_parts = var["ALT"].split(":")
                alt = (alt_parts[2] if len(alt_parts) > 2 else var["ALT"]).replace(">", "")
                event = str(var["INFO"].get("INTERNAL", "")).replace(",", "_")
                meinfo = str(var["INFO"].get("MEINFO", "")).replace(",", " ")
                custom = [
                    f"SCOUT_CUSTOM=Repeat Element|{alt}",
                    f"Subtype|{meinfo}",
                    f"Consequence|{event}",
                    f"Target Site Duplication|{var['INFO'].get('TSD', '')}",
                    f"Filter|{qcval.get(var['FILTER'], '')}",
                    f"Assess|{assess.get(str(var['INFO'].get('ASSESS', '')), '')}",
                ]
                gt = var["GT"][0]
                depth = leading_float(gt.get("DP"))
                alt_depth = leading_float(gt.get("AD"))
                vaf_value = f"{alt_depth / depth:.4f}" if depth else "0.0000"
                row = [
                    var["CHROM"],
                    var["POS"],
                    var["ID"],
                    var["REF"],
                    var["ALT"],
                    var["QUAL"],
                    "PASS",
                    f"SVTYPE=INS;END={var['POS']};{','.join(custom)}",
                    "GT:VAF:VD:DP",
                    f"{gt.get('GT', '')}:{vaf_value}:{gt.get('AD', '')}:{gt.get('DP', '')}",
                ]
                print("\t".join(row), file=out_fh)

    out_path = out_file or f"{sample_id}.melt.merged.vcf"
    if not mod_files:
        with open(vcf_header) as header_fh, open(out_path, "w") as out_fh:
            for line in header_fh:
                out_fh.write(line.replace("FORMAT", f"FORMAT\t{sample_id}") if line.startswith("#CHROM") else line)
        return out_path
    if shutil.which("vcf-concat") and shutil.which("vcf-sort"):
        with open(out_path, "w") as out_fh:
            concat = subprocess.Popen(["vcf-concat", *mod_files], stdout=subprocess.PIPE, text=True)
            subprocess.run(["vcf-sort", "-c"], stdin=concat.stdout, stdout=out_fh, text=True, check=True)
            concat.wait()
    else:
        lines = []
        header = []
        for idx, path in enumerate(mod_files):
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        if idx == 0:
                            header.append(line)
                    else:
                        lines.append(line)
        lines.sort(key=lambda row: (_natural_key(row.split("\t", 1)[0]), int(row.split("\t", 2)[1])))
        with open(out_path, "w") as out_fh:
            out_fh.writelines(header + lines)
    return out_path


def overlapping_genes(in_bed, genes_bed, out_file):
    """Annotate BED intervals with overlapping gene names using bedtools."""
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        tmp_path = tmp.name
        with open(in_bed) as in_fh:
            for line in in_fh:
                if line.startswith("@") or line.startswith("CONTIG") or line.startswith("REF"):
                    continue
                tmp.write(line)
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", tmp_path, "-b", genes_bed, "-loj"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
    finally:
        Path(tmp_path).unlink(missing_ok=True)
    with open_output(out_file) as out_fh:
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            out_fh.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[-1]}\n")


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


def coyote_segmentator(vcf_file, sample_id, panel_files, genes_bed, normal=False, out_prefix=None):
    """Create coyote CN segment BED files from a merged SV VCF."""
    panel_genes = {}
    for panel in panel_files.split(","):
        _read_panel(panel, panel_genes)

    out_prefix = out_prefix or sample_id
    records = []
    with vcf.VCFReader(vcf_file) as reader:
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
            records.append(f"{var['CHROM']}\t{start}\t{end}\t{probes}\t{fold}\t{sign}\t{callers}\t{pr}:{sr}")

    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        tmp_path = tmp.name
        for row in records:
            tmp.write(row + "\n")
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", tmp_path, "-b", genes_bed, "-loj"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    filtered_path = f"{out_prefix}.cn-segments.panel.bed"
    unfiltered_path = f"{out_prefix}.cn-segments.bed"
    groups = []
    previous_region = None
    genes = []
    previous_bed = None
    for line in result.stdout.splitlines():
        fields = line.split("\t")
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


def build_filter_gatk_parser():
    """Build parser for GATK SV conversion."""
    parser = ArgumentParser(description="Convert GATK CNV records to SV-style VCF records.")
    parser.add_argument("--vcf", required=True, help="Input GATK VCF.")
    parser.add_argument("--out", required=True, help="Output VCF.")
    return parser


def main_filter_gatk(argv=None):
    """Run GATK SV conversion."""
    args = build_filter_gatk_parser().parse_args(argv)
    filter_gatk(args.vcf, args.out)
    return 0


def build_parser():
    parser = ArgumentParser(description="Annotate BED intervals with overlapping genes.")
    parser.add_argument("--bed", required=True, help="Input BED file.")
    parser.add_argument("--genes-bed", required=True, help="Gene BED file.")
    parser.add_argument("--out", required=True, help="Output table.")
    return parser

def main(argv=None):
    args = build_parser().parse_args(argv)
    overlapping_genes(args.bed, args.genes_bed, args.out)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
