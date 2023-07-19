#!/usr/bin/env nextflow

include { PREPROCESSINTERVALS      } from '../../modules/local/GATK_references/main'
include { COUNT_READS              } from '../../modules/local/GATK_references/main'
include { ANNOTATE_GC              } from '../../modules/local/GATK_references/main'
include { CORRECT_GC               } from '../../modules/local/GATK_references/main'
include { REMOVE_ALT_CONTIGS       } from '../../modules/local/GATK_references/main'
include { SCATTER_INTERVALS        } from '../../modules/local/GATK_references/main'

workflow BED_INTERVALS {
	take:
		prefix         // val(prefix) name of reference
		sample		   // val(id), file(cram), file(crai), file(bai)

	main:
        PREPROCESSINTERVALS(prefix)
		ANNOTATE_GC(PREPROCESSINTERVALS.out.preprocessed)
		COUNT_READS(sample,PREPROCESSINTERVALS.out.preprocessed)
		CORRECT_GC(PREPROCESSINTERVALS.out.preprocessed.join(ANNOTATE_GC.out.annotated_intervals), COUNT_READS.out.count_tsv.groupTuple())
		REMOVE_ALT_CONTIGS(CORRECT_GC.out.corrected_intervals)
		if (params.panel) {
			scattered = REMOVE_ALT_CONTIGS.out.noaltcontigs
		}
		else {
			SCATTER_INTERVALS(REMOVE_ALT_CONTIGS.out.noaltcontigs, params.scatter_size)  // scatter bins. 380000 per bin creates 7 scatters for human wgs, this is used to fasten up calling downstream
			scattered = SCATTER_INTERVALS.out.scatter
		}
		
		

	emit:
        intervals = REMOVE_ALT_CONTIGS.out.noaltcontigs
		intervals_scattered = scattered
		counts = COUNT_READS.out.count_tsv

}