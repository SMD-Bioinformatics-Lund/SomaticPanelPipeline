#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
    take:
        csv     // file(csv)

    main:
        CSV_CHECK(csv)
        checkedCSV = CSV_CHECK.out.csv.splitCsv(header:true, sep:',').set { csvmap }

        reads = csvmap.map { row -> create_fastq_channel(row) }
            .collect(flat: false)
            .flatMap { samples ->
                def ids = samples.collect { it[1].id }
                if (ids.unique(false).size() != ids.size()) {
                    error "Sample sheet contains duplicate sample IDs"
                }
                def groups = samples.groupBy { it[0] }
                groups.each { group, members ->
                    def types = members.collect { it[1].type?.toString()?.toLowerCase() }
                    if (!group || types.any { !(it in ['t', 'tumor', 'n', 'normal']) }) {
                        error "Each sample must have a group and a tumor/normal type"
                    }
                    if (members.size() > 2 || (members.size() == 2 &&
                        !(types.count { it in ['t', 'tumor'] } == 1 && types.count { it in ['n', 'normal'] } == 1))) {
                        error "Group ${group} must contain one sample or one tumor and one normal"
                    }
                }
                samples.collect { sample ->
                    def meta = sample[1] + [paired: groups[sample[0]].size() == 2]
                    tuple(sample[0], meta, sample[2], sample[3])
                }
            }

        // FASTQ
        fastq = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            (read1.endsWith("fastq.gz") || read1.endsWith("fq.gz")) &&
            (read2.endsWith("fastq.gz") || read2.endsWith("fq.gz"))
        }

        // BAM + BAI
        bam = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            read1.endsWith("bam") && (read2.endsWith("bai") || read2.endsWith("bam.bai"))
        }

        // VCF + index
        vcf = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            read1.endsWith("vcf") &&
                (read2.endsWith("tbi") || read2.endsWith("csi") || read2.endsWith("vcf.gz.tbi"))
        }

        meta = reads
            .map { group, meta, r1, r2 ->
                tuple( group, meta)
            }

    emit:
        fastq
        bam
        vcf
        meta
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
	// create meta map
	def meta = [:]
	meta.id             	= row.id
	meta.group              = row.group
	meta.diagnosis          = row.diagnosis
	meta.type               = row.type
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = row.containsKey("ffpe") && row.ffpe ? true : false
	meta.purity_raw         = (row.containsKey("purity") ? row.purity : false)
	meta.purity             = normalize_purity(meta.purity_raw)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
    meta.sex                = (row.containsKey("sex") ? row.sex : false)
	meta.clarity_pool_id    = row.clarity_pool_id
    def sub = false
    if (meta.reads && params.sample) {
        if (meta.reads.toInteger() > params.sample_val) {
            sub = (params.sample_val / meta.reads.toInteger()).round(2)
            if (sub == 1.00) sub = 0.99
        }
    }
    meta.sub = sub
	// add path(s) of the fastq file(s) to the meta map
	def fastq_meta = []
	fastq_meta = [row.group, meta, file(row.read1), file(row.read2) ]

	return fastq_meta
}

def normalize_purity(purity) {
    def value = purity?.toString()?.trim()
    if (!value || value ==~ /(?i:false|0+(?:\.0*)?)/) {
        return false
    }
    return purity
}
