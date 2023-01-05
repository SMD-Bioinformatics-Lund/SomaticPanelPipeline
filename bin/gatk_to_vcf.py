import argparse
from datetime import date


def print_header(sample):
    date2d = str(date.today())
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write("##fileDate=%s\n" % date2d)
    vcf_out.write("##source=GATK v4.1.0\n")
    vcf_out.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    vcf_out.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
    vcf_out.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_out.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_out.write("##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description=\"Fold change\">\n")
    vcf_out.write("##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description=\"Log fold change\">\n")
    vcf_out.write("##INFO=<ID=PROBES,Number=1,Type=Integer,Description=\"Number of probes in CNV\">\n")
    vcf_out.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    vcf_out.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    vcf_out.write("##ALT=<ID=CNV,Description=\"Copy number variable region\">\n")
    vcf_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_out.write("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n")
    vcf_out.write("##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description=\"Copy number genotype quality for imprecise events\">\n")
    vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)




parser = argparse.ArgumentParser()
parser.add_argument('-s', '--segments')
parser.add_argument('-o', '--out')
parser.add_argument('-p', '--tumorpurity')
args = parser.parse_args()
seg_in = open(args.segments)

vcf_out = open(str(args.out), "w")

for line in seg_in:
    line = line.strip().split("\t")
    
    
    if line[0] == "@RG":
        sample = line[2].split(":")
        sample = sample[1]
    elif line[0][0] == "@":
        continue
    elif line[0] == "CONTIG":
        print_header(sample)
    else:
        variant = []
        chrom = line[0]
        variant.append(chrom)
        start = line[1]
        variant.append(start)
        end = line[2]
        num_probes = line[3]
        log2 = float(line[4])
        gatk_call = line[5]

        # ID #
        id = "CNV_" + str(chrom) + "_" + str(start) + "_" + str(end)
        variant.append(id)
 
        ## Copy Number, and adjusted for tumor purity
        purity = args.tumorpurity
        if purity == "false":
            purity = 1
        fold_change = round(pow(2, (float(log2)) ),6)
        if args.tumorpurity:
            cn = round(2 * fold_change * 1/float(purity)) 
        else:
            cn = round(2 * fold_change)

        variant.append("N")
        # SVTYPE #
        if gatk_call == "+":
            svtype = "DUP"
            variant.append("<%s>" % svtype)
        elif gatk_call == "-":
            svtype = "DEL"
            variant.append("<%s>" % svtype)
        else:
            svtype = "."
            continue
        
        ## INFO-FIELD
        info = []
        ## length
        svlen = int(end) - int(start) + 1
        if svtype == "DEL":
            svlen = 0 - svlen
        
        info.append("IMPRECISE")
        info.append("SVLEN=%s" % svlen)
        info.append("SVTYPE=%s" % svtype)
        info.append("END=%s" % end)
        info.append("FOLD_CHANGE_LOG=%s" % log2)
        info.append("FOLD_CHANGE=%s" % fold_change)
        info.append("PROBES=%s" % num_probes)
        
        variant.append(".\t.") #filter and qual
        variant.append(';'.join(info))
        
        ## FORMAT and sample FIELD
        format = []
        gt = "0/1"
        if cn < 0.2:
            gt = "1/1"
        format.append("GT:CN:CNQ")
        format.append("%s:%s:%s" % (gt, cn, num_probes))

        variant.append('\t'.join(format))


        vcf_out.write("\t".join(variant))
        vcf_out.write("\n")
        #print("\t".join(variant))
