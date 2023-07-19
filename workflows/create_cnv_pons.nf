#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { CHECK_INPUT                     } from '../subworkflows/local/create_input_cnvrefs'
include { BED_INTERVALS                   } from '../subworkflows/local/bed_intervals'
include { CALL_COHORT                     } from '../subworkflows/local/call_cohort'
include { CNVKITREFS                      } from '../subworkflows/local/cnvkit_refs'

println(params.genome_file)

csv = file(params.csv)
println(csv)

workflow CREATE_REF {

	// Checks input, creates meta-channel and decides whether data should be downsampled //
	CHECK_INPUT ( Channel.fromPath(csv) )
    CNVKITREFS( 
        CHECK_INPUT.out.sample,
        params.name
    )
	BED_INTERVALS(
		params.name,
		CHECK_INPUT.out.sample
	)
	.set { ch_bed }
	CALL_COHORT (
		ch_bed.intervals,
		ch_bed.intervals_scattered,
		ch_bed.counts
	)
         
		 
		 
}

workflow {
	CREATE_REF()
}

