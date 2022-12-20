#!/usr/bin/env nextflow

include { SEQTK                } from '../../modules/local/seqtk/main'
include { FASTP                } from '../../modules/local/fastp/main'


workflow SAMPLE {
	take:
		fastq

	main:
		// Only sub-sample if meta.sub == value
		SEQTK ( fastq.filter{ item -> item[1].sub != false } )
		//combine with any sample that does not get sub-sampled
		fastq_sample = SEQTK.out.fastq_sub.mix( fastq.filter{ item -> item[1].sub == false } )

		if (params.trimfq) {
			FASTP ( fastq_sample )
			fastq_done = FASTP.out.fastq_trimmed
		}
		else {
			fastq_done = fastq_sample
		}

	emit:
		fastq_trim = fastq_done

}