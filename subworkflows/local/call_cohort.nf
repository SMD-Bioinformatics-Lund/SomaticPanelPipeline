#!/usr/bin/env nextflow

include { COHORT_PLOIDY           } from '../../modules/local/GATK_references/main'
include { COHORT_CALL             } from '../../modules/local/GATK_references/main'
include { COHORT_CALL_PANEL       } from '../../modules/local/GATK_references/main'
include { GATK_SOM_PON            } from '../../modules/local/GATK_references/main'


workflow CALL_COHORT {
	take:
		intervals           // val(prefix), file(intervals) name of reference
        intervals_scattered // val(prefix), file(scatters.tar)
		counts              // val(prefix), file(id), file(tsv)

	main:
        COHORT_PLOIDY(intervals,counts.groupTuple())
		scatters = intervals_scattered.map{ val-> val[1] }.flatten()
		if (params.panel) {
			COHORT_CALL_PANEL(counts.groupTuple().join(COHORT_PLOIDY.out.ploidy).join(intervals))
		}
		else {
			COHORT_CALL(counts.groupTuple().join(COHORT_PLOIDY.out.ploidy).combine(scatters))
		}
		GATK_SOM_PON( counts.groupTuple() )
		
	emit:
        ploidy = COHORT_PLOIDY.out.ploidy


}