#!/usr/bin/env python3
"""Sample fingerprinting from BAM pileups."""

import glob
from argparse import ArgumentParser
from pathlib import Path

DEFAULT_SNPBED = Path(__file__).resolve().parent / "HPA_1000G_final_38.bed"
DEFAULT_XYBED = Path(__file__).resolve().parent / "xy_38.bed"


def read_bed(path):
    """Read SNP marker BED data."""
    variants = {}
    chromosomes = {}
    allele_freq = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            loc = f"{cols[0]}:{cols[2]}"
            chromosomes[cols[0]] = chromosomes.get(cols[0], 0) + 1
            variants[loc] = cols[3] if len(cols) > 3 and cols[3] else loc
            if len(cols) > 4 and cols[4]:
                allele_freq[loc] = cols[4]
    return variants, chromosomes, allele_freq


def read_sample_metadata(path, nocheck=False):
    """Read BAM metadata table."""
    files = {}
    used_names = set()
    unknown_count = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if not Path(cols[0]).is_file() and not nocheck:
                raise FileNotFoundError(f"ERROR: File not found '{cols[0]}'")
            unknown_count += 1
            name = cols[1] if len(cols) > 1 and cols[1] else f"unknown{unknown_count}"
            if name in used_names:
                raise ValueError(f"ERROR: Sample IDs must be unique ({name} found multiple times)")
            used_names.add(name)
            files[cols[0]] = {
                "name": name,
                "individual": cols[2] if len(cols) > 2 else "",
                "sex": cols[3] if len(cols) > 3 else "",
            }
    return files


def get_bampaths(mask, nocheck=False):
    """Resolve BAM paths from a metadata file or glob mask."""
    files = sorted(glob.glob(mask))
    if not files:
        raise FileNotFoundError(f"No file(s) found: {mask}")
    if len(files) == 1:
        path = files[0]
        if not Path(path).is_file() and not nocheck:
            raise FileNotFoundError(f"No file(s) found: {mask}")
        if path.endswith((".bam", ".sam")):
            return {path: {"name": path, "individual": "", "sex": ""}}
        return read_sample_metadata(path, nocheck=nocheck)
    return {path: {"name": path, "individual": "", "sex": ""} for path in files}


def count_bases(base_text):
    """Count A/C/G/T observations in a mpileup base string."""
    if not base_text:
        return {"A": 0, "C": 0, "G": 0, "T": 0}
    upper = base_text.upper()
    return {base: upper.count(base) for base in "ACGT"}


def merge_files(files, out_path):
    """Concatenate files into one output file."""
    with open(out_path, "w") as out_fh:
        for path in files:
            if not path:
                continue
            with open(path) as in_fh:
                out_fh.writelines(in_fh)


def get_base_freqs_from_bams(file_data, pileup_file):
    """Collect base counts by locus and sample from a precomputed mpileup file."""
    bams = sorted(file_data)
    freqs = {}
    with open(pileup_file) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            loc = f"{cols[0]}:{cols[1]}"
            freqs.setdefault(loc, {"samples": {}})
            idx = 3
            for bam in bams:
                sample_name = file_data[bam]["name"]
                depth_text, bases_text, _qual = cols[idx : idx + 3]
                bases = count_bases(bases_text)
                freqs[loc]["samples"][sample_name] = {
                    "bases": bases,
                    "depth": sum(bases.values()) or int(depth_text),
                }
                idx += 3
    return freqs


def largest_base(counts, position):
    """Return the base with the nth largest count."""
    ordered = sorted(counts, key=lambda base: counts[base], reverse=True)
    return ordered[position - 1]


def do_genotyping(data, loci, allele_freq, observed_af=False):
    """Assign genotypes to each sample at each marker."""
    for loc in loci:
        if loc not in data:
            data[loc] = {"samples": {}, "Ntotal": 0, "Ncallable": 0}
            continue
        genotype_counts = {}
        callable_count = 0
        total_count = 0
        for sample_id in sorted(data[loc]["samples"]):
            sample = data[loc]["samples"][sample_id]
            counts = sample["bases"]
            total = sample["depth"]
            max_base = largest_base(counts, 1)
            second_base = largest_base(counts, 2)
            if total >= 6:
                allele_proportion = (counts[max_base] - counts[second_base]) / total
                if allele_proportion < 0.6:
                    genotype = "/".join(sorted((max_base, second_base)))
                elif allele_proportion > 0.9:
                    genotype = max_base
                else:
                    genotype = "unclear"
            else:
                genotype = "lowdata"
            if genotype not in {"unclear", "lowdata"}:
                sample["basecall"] = genotype
                genotype_counts[genotype] = genotype_counts.get(genotype, 0) + 1
                callable_count += 1
            else:
                sample["uncalled"] = genotype
            total_count += 1
        data[loc]["Ntotal"] = total_count
        data[loc]["Ncallable"] = callable_count
        if not callable_count:
            data[loc]["skip"] = True
            continue
        if allele_freq.get(loc) and not observed_af:
            alt, alt_af = allele_freq[loc].split("/", 1)
            data[loc]["alt"] = alt
            data[loc]["alt_af"] = alt_af


def mean(values):
    """Return the arithmetic mean."""
    return sum(values) / len(values) if values else 0


def determine_sex(data, sex_loci):
    """Predict sex using XY marker coverage."""
    if not sex_loci:
        return {}
    avg_cov = {}
    for loc, item in data.items():
        if loc in sex_loci:
            continue
        for sample_id, sample in item["samples"].items():
            avg_cov.setdefault(sample_id, []).append(sample["depth"])
    avg_cov = {sample_id: mean(values) for sample_id, values in avg_cov.items()}
    sample_data = {}
    male_hits = {}
    female_hits = {}
    for loc in list(sex_loci):
        if loc not in data:
            continue
        for sample_id, sample in data[loc]["samples"].items():
            ratio = sample["depth"] / avg_cov.get(sample_id, 1)
            if ratio >= 0.2:
                male_hits[sample_id] = male_hits.get(sample_id, 0) + 1
            elif ratio <= 0.05:
                female_hits[sample_id] = female_hits.get(sample_id, 0) + 1
        del data[loc]
    for sample_id in set(male_hits) | set(female_hits):
        if male_hits.get(sample_id) and not female_hits.get(sample_id):
            pred = "M"
        elif female_hits.get(sample_id) and not male_hits.get(sample_id):
            pred = "F"
        else:
            pred = "unclear"
        sample_data.setdefault(sample_id, {"sex": {}})["sex"]["pred"] = pred
    return sample_data


def plink_gt(basecall):
    """Convert a basecall to PLINK two-column genotype format."""
    if not basecall or basecall in {"unclear", "lowdata"}:
        return "0\t0"
    if len(basecall) == 1:
        return f"{basecall}\t{basecall}"
    return "\t".join(basecall.split("/", 1))


def print_genotype_table(data, sample_data, out_prefix, annotation, variants, long=False, position=False):
    """Write genotype, statistics, and sex prediction tables."""
    samples = sorted({sample for loc in data.values() for sample in loc.get("samples", {})})
    with open(f"{out_prefix}.genotypes", "w") as gt_fh:
        if not long:
            gt_fh.write("loc" + "".join(f"\t{sample}" for sample in samples) + "\n")
        total_callable = total_sites = 0
        for loc in sorted(data):
            snp_id = loc if position else variants.get(loc, loc)
            if not long:
                gt_fh.write(snp_id)
            for sample_id in samples:
                genotype = data[loc].get("samples", {}).get(sample_id, {}).get("basecall")
                if long:
                    gt_fh.write(f"{sample_id}\t{sample_id}\t{snp_id}\t{plink_gt(genotype)}\n")
                else:
                    gt_fh.write(f"\t{genotype if genotype is not None else 'NA'}")
            if not long:
                gt_fh.write("\n")
            total_callable += data[loc].get("Ncallable", 0)
            total_sites += data[loc].get("Ntotal", 0)
    with open(f"{out_prefix}.stats", "w") as stats_fh:
        stats_fh.write(
            f"Average callability: {100 * (total_callable / total_sites):.2f}%\n"
            if total_sites
            else "Average callability: 0.00%\n"
        )
    with open(f"{out_prefix}.sex", "w") as sex_fh:
        sex_fh.write("sample\tpredicted_sex\tannotated_sex\n")
        for bam in sorted(annotation):
            sample_id = annotation[bam]["name"]
            pred = sample_data.get(sample_id, {}).get("sex", {}).get("pred", "-")
            sex_fh.write(f"{sample_id}\t{pred}\t{annotation[bam].get('sex') or '-'}\n")


def fingerprint_bams(
    bam,
    out_prefix,
    pileup_file=None,
    bed=str(DEFAULT_SNPBED),
    bedxy=str(DEFAULT_XYBED),
    nocheck=False,
    long=False,
    position=False,
    obsaf=False,
):
    """Create sample fingerprint output files from a precomputed mpileup."""
    meta_data = get_bampaths(bam, nocheck=nocheck)
    variants, chromosomes, allele_freq = read_bed(bed)
    sex_loci = {}
    if bedxy and bedxy != "none":
        sex_loci, xy_chromosomes, _xy_af = read_bed(bedxy)
        chromosomes.update(xy_chromosomes)
    pileup_file = pileup_file or f"{out_prefix}.pile.tmp"
    data = get_base_freqs_from_bams(meta_data, pileup_file)
    sample_data = determine_sex(data, sex_loci)
    do_genotyping(data, variants, allele_freq, observed_af=obsaf)
    print_genotype_table(data, sample_data, out_prefix, meta_data, variants, long=long, position=position)


def build_parser():
    """Build the provider parser."""
    parser = ArgumentParser(description="Create sample fingerprints from BAM pileups.")
    parser.add_argument("--bam", required=True, help="BAM file, metadata table, or glob mask.")
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--pileup", help="Precomputed samtools mpileup output. Defaults to <out>.pile.tmp.")
    parser.add_argument("--bed", default=str(DEFAULT_SNPBED), help="SNP marker BED.")
    parser.add_argument("--bedxy", default=str(DEFAULT_XYBED), help="Sex marker BED, or 'none'.")
    parser.add_argument("--nocheck", action="store_true", help="Skip file/chromosome checks.")
    parser.add_argument("--long", action="store_true", help="Write long-format genotype table.")
    parser.add_argument("--position", action="store_true", help="Use genomic positions as marker IDs.")
    parser.add_argument("--obsaf", action="store_true", help="Use observed allele frequencies.")
    return parser


def main(argv=None):
    """Run sample fingerprinting."""
    args = build_parser().parse_args(argv)
    fingerprint_bams(
        args.bam, args.out, args.pileup, args.bed, args.bedxy, args.nocheck, args.long, args.position, args.obsaf
    )
    return 0
