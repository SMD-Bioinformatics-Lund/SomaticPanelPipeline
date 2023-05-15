import json
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-j', '--json')
parser.add_argument('-o', '--out')
parser.add_argument('-i', '--id')
args = parser.parse_args()
seg_in = open(args.json)


genefuse = open(args.json)
id = args.id
vcf_out = open(str(args.out), "w")


vcf_out.write("##fileformat=VCFv4.1\n")
vcf_out.write("##source=genefuse\n")
vcf_out.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
vcf_out.write("##INFO=<ID=UR,Number=1,Type=Integer,Description=\"Number of unique reads supporting variant\">\n")
vcf_out.write("##FORMAT=<ID=UR,Number=.,Type=Integer,Description=\"Unique reads supporting fusion\">\n")


reverse = {
    'T': 'A',
    'A': 'T',
    'G' : 'C',
    'C': 'G'
}
gf_dict = json.load(genefuse)
command = gf_dict['command']
command  = "##cmdline=" + command
vcf_out.write("%s\n" % command)
vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % id)


for fusion in gf_dict['fusions']:

    # left info
    chrom = gf_dict['fusions'][fusion]['left']['gene_chr']
    pos = gf_dict['fusions'][fusion]['left']['position']
    gene_name = gf_dict['fusions'][fusion]['left']['gene_name']
    strand = gf_dict['fusions'][fusion]['left']['strand']
    reference = gf_dict['fusions'][fusion]['left']['reference'][-1]
    
    if strand == "reversed":
        reference = reverse[reference]
    neg = 0
    if pos < 0:
        neg = 1
        pos = pos * -1
    # right info
    r_chrom = gf_dict['fusions'][fusion]['right']['gene_chr']
    r_pos = gf_dict['fusions'][fusion]['right']['position']
    r_gene_name = gf_dict['fusions'][fusion]['right']['gene_name']
    if r_pos < 0:
        r_neg = 1
        r_pos = r_pos * -1
    # off-by-one-error-correction
    r_pos = r_pos + 1
    if neg:
        alt = "[" + str(r_chrom) + ":" + str(r_pos) + "[" + reference
    else:
        alt = reference + "]" + str(r_chrom) + ":" + str(r_pos) + "]"
    # off-by-one-error-correction
    pos = pos + 1
    vcf_out.write("%s\t" % chrom)
    vcf_out.write("%s\t" % pos)
    vcf_out.write("genefuse\t")
    vcf_out.write("%s\t" % reference)
    vcf_out.write("%s\t" % alt)
    vcf_out.write(".\t")
    vcf_out.write("PASS\t")
    vcf_out.write("SVTYPE=BND;")
    ur = "UR="+str(gf_dict['fusions'][fusion]['unique'])
    vcf_out.write("%s\t" % ur)
    vcf_out.write("UR\t")
    vcf_out.write("%s\n" % str(gf_dict['fusions'][fusion]['unique']))
    