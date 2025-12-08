#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
    take:
        csv     // file(csv)
        paired  // boolean

    main:
        CSV_CHECK(csv)
        checkedCSV = CSV_CHECK.out.csv.splitCsv(header:true, sep:',').set { csvmap }

        reads     = csvmap.map { create_fastq_channel(it, paired) }

        // FASTQ
        ch_fastq = reads.filter { it ->
            def r1 = it[2].toString()
            def r2 = it[3].toString()
            r1.endsWith("fastq.gz") || r1.endsWith("fq.gz") &&
            r2.endsWith("fastq.gz") || r2.endsWith("fq.gz")
        }

        // BAM + BAI
        ch_bam = reads.filter { it ->
            def r1 = it[2].toString()
            def r2 = it[3].toString()
            r1.endsWith("bam") && (r2.endsWith("bai") || r2.endsWith("bam.bai"))
        }

        // VCF + index
        ch_vcf = reads.filter { it ->
            def r1 = it[2].toString().toLowerCase()
            def r2 = it[3].toString().toLowerCase()
            r1.endsWith("vcf") &&
                (r2.endsWith("tbi") || r2.endsWith("csi") || r2.endsWith("vcf.gz.tbi"))
        }

        meta = csvmap.map { create_samples_channel(it, paired) }

    emit:
        ch_fastq
        ch_bam
        ch_vcf
        meta
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, paired) {
	// create meta map
	def meta = [:]
	meta.id             	= row.id
	meta.group              = row.group
	meta.diagnosis          = row.diagnosis
	meta.type               = row.type
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = row.containsKey("ffpe") && row.ffpe ? true : false
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
	meta.clarity_pool_id    = row.clarity_pool_id
    meta.paired             = paired
	sub = false
	if (meta.reads && params.sample) {  
        if (meta.reads.toInteger() > params.sample_val) {
        	sub = (params.sample_val / meta.reads.toInteger()).round(2)
        	if (sub == 1.00){
                sub = 0.99
        	}
        }
        else {
        	sub = false
        }
	}
	meta.sub = sub
	// add path(s) of the fastq file(s) to the meta map
	def fastq_meta = []
	fastq_meta = [row.group, meta, file(row.read1), file(row.read2) ]

	return fastq_meta
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row, paired) {
	def meta                = [:]
	meta.id                 = row.id
	meta.group              = row.group
	meta.diagnosis          = row.diagnosis
	meta.type               = row.type
	meta.clarity_sample_id  = row.clarity_sample_id
	meta.ffpe               = row.containsKey("ffpe") && row.ffpe ? true : false
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	meta.sequencing_run     = row.sequencing_run
	meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
	meta.clarity_pool_id    = row.clarity_pool_id
    meta.paired             = paired

	def sample_meta = [row.group, meta]
	return sample_meta
}
