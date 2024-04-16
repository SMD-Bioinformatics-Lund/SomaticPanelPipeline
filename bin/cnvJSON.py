#!/usr/bin/env python
import logging
import json
import sys
from pprint import pprint
"""
Each variant needs this format. Could have less information depending on variantcaller
manta-like specific
PR and SR 

read-depth based specific
nprobes
ratio

non-mandatory
callers

{
        "_id" : ObjectId("6537d24d9bc472e74b7df46e"),
        "SR" : 0,
        "PR" : 0,
        "chr" : 7,
        "nprobes" : 87,
        "genes" : [
                {
                        "gene" : "ADCK2"
                },
                {
                        "gene" : "NDUFB2"
                },
                {
                        "gene" : "DENND2A"
                },
                {
                        "gene" : "BRAF",
                        "cnv_type" : "unspec",
                        "class" : "somatic"
                }
        ],
        "size" : 344724,
        "callers" : "cnvkit",
        "SAMPLE_ID" : "6537d23c9bc472e74b7de061",
        "end" : 140924045,
        "ratio" : 2.95721,
        "start" : 140579321
}
"""

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

def read_bedlike_vcf(file_name,panel_genes):
    variants = {}
    try:
        with open(file_name, 'r') as vcf:
            for line in vcf:
                line = line.split("\t")
                chrom = line[0]
                start = line[1]
                genematch = line.pop().rstrip('\n')
                LOG.debug("Matched gene: %s", genematch)
                if not genematch:
                    genematch = ["."]
                info = info_field(line[7])
                gt = gt_field(line[8],line[9])
                end = info["END"]
                varid = chrom+":"+start+"-"+end
                LOG.debug("Processing varid: %s", varid)
                varinfo = get_varinfo(info,gt)
                varinfo["chr"] = str(chrom)
                varinfo["start"] = int(start)
                varinfo["end"] = int(end)

                if varid in variants:
                    genes = [genematch]
                    for gene in variants[varid]["genes"]:
                        genes.append(gene["gene"])
                    LOG.debug("Sending list of genes: %s", genes)
                    gene_list,has_match = collate_genes(genes,panel_genes)
                    variants[varid]["genes"] = gene_list
                    if not "has_match" in variants[varid] and has_match:
                        variants[varid]["has_match"] = True
                else:
                    genes = [genematch]
                    genes_list,has_match = collate_genes(genes,panel_genes)
                    LOG.debug("Creating varinfo['genes'] from genes_list: %s", genes_list)
                    varinfo["genes"] = genes_list
                    variants[varid] = varinfo
                    if not "has_match" in variants[varid] and has_match:
                        variants[varid]["has_match"] = True
            
    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found.")
    return variants
    
def info_field(info):
    infodict = {}
    info = info.split(";")
    for item in info:
        try:
            key,value = item.split('=')
            infodict[key] = value
        except:
            key = item
            value = True
            infodict[key] = value
    return infodict

def get_varinfo(info,gt):
    varinfo = {}
    if "PROBES" in info:
        varinfo["nprobes"] = int(info["PROBES"])
    if "set" in info:
        varinfo["callers"] = info["set"]
    if "FOLD_CHANGE_LOG" in info:
        varinfo["ratio"] = float(info["FOLD_CHANGE_LOG"])
    if "SVLEN" in info:
        varinfo["size"] = int(info["SVLEN"])
    if "PR" in gt:
        varinfo["PR"] = gt["PR"]
    if "SR" in gt:
        varinfo["SR"] = gt["SR"]
    return varinfo

def gt_field(gt_keys,gt_values):
    gt_dict = {}
    gt_keys = gt_keys.split(':')
    gt_keys = [item.rstrip('\n') for item in gt_keys]
    gt_values = gt_values.split(':')
    gt_values = [item.rstrip('\n') for item in gt_values]
    for i in range(len(gt_keys)):
        gt_dict[gt_keys[i]] = gt_values[i]
    return gt_dict

def read_panel_genes(panel_file):
    panel_genes = {}
    with open(panel_file, 'r') as panel:
        for line in panel:
            line = line.split("\t")
            tmp = {}
            tmp["cnv_type"] = line[2].rstrip('\n')
            tmp["class"] = line[1]
            panel_genes[line[0]] = tmp
    return panel_genes


def collate_genes(genes,panel_genes):
    genes_list = []
    has_match = False
    LOG.debug("Processing genes: %s", genes)
    for gene in genes:
        if gene in panel_genes:
            tmp = {}
            tmp["gene"] = gene
            tmp["class"] = panel_genes[gene]["class"]
            tmp["cnv_type"] = panel_genes[gene]["cnv_type"]
            genes_list.append(tmp)
            LOG.debug("Panel gene match data: %s", genes_list)
            has_match = True
        else:
            tmp = {}
            tmp["gene"] = gene
            genes_list.append(tmp)
    return genes_list,has_match

file_name = sys.argv[1]
panel_file = sys.argv[2]
sample_id = sys.argv[3]
panel_genes = read_panel_genes(panel_file)
variants = read_bedlike_vcf(file_name,panel_genes)
json_all = sample_id+"_cnvs.json"
json_matched = sample_id+"_cnvs_panelmatched.json"
variants_matched = {}
for var in variants:
    if "has_match" in variants[var]:
        del variants[var]["has_match"]
        variants_matched[var] = variants[var]
with open(json_all, 'w') as json_file:
    json.dump(variants, json_file)

with open(json_matched, 'w') as json_file:
    json.dump(variants_matched, json_file)