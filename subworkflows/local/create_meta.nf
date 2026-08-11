#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK } from '../../modules/local/check_input/main'
include { UNSPRING  } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
    take:
        csv     // file(csv)
        paired  // boolean

    main:
        CSV_CHECK(csv)
        checkedCSV = CSV_CHECK.out.csv.splitCsv(header:true, sep:',').set { csvmap }

        reads     = csvmap.map { row -> create_fastq_channel(row, paired) }

        // FASTQ
        fastq_input = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            if (!read1 || !read2) return false
            (read1.endsWith("fastq.gz") || read1.endsWith("fq.gz")) &&
            (read2.endsWith("fastq.gz") || read2.endsWith("fq.gz"))
        }

        // SPRING archive
        spring = reads
            .filter { group, meta, r1, r2 ->
                def read1 = r1.toString()
                def read2 = r2.toString()
                read1.endsWith(".spring") || read2.endsWith(".spring")
            }
            .map { group, meta, r1, r2 ->
                def spring = r1.toString().endsWith(".spring") ? r1 : r2
                tuple( group, meta, spring )
            }

        UNSPRING ( spring )

        fastq = fastq_input.mix( UNSPRING.out.fastq )

        // BAM + BAI
        bam = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            if (!read1 || !read2) return false
            read1.endsWith("bam") && (read2.endsWith("bai") || read2.endsWith("bam.bai"))
        }

        // VCF + index
        vcf = reads.filter { group, meta, r1, r2 ->
            def read1 = r1.toString()
            def read2 = r2.toString()
            if (!read1 || !read2) return false
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
    meta.sex                = (row.containsKey("sex") ? row.sex : false)
	meta.clarity_pool_id    = row.clarity_pool_id
    meta.paired             = paired
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
	fastq_meta = [row.group, meta, create_input_file(row.read1), create_input_file(row.read2) ]

	return fastq_meta
}

def create_input_file(value) {
    def path = value?.toString()?.trim()
    return path ? file(path) : ''
}
