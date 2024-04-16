#!/usr/bin/env nextflow

include { ANNOTATE_VEP                          } from '../../modules/local/filters/main'
include { COYOTE_SEGMENTS                       } from '../../modules/local/filters/main'
include { MERGE_SEGMENTS                        } from '../../modules/local/filters/main'
include { INTERSECT as GENE_INTERSECT           } from '../../modules/local/bedtools/main'
include { COYOTE_SEGMENTS_JSON                  } from '../../modules/local/filters/main'

workflow CNV_ANNOTATE {
	take: 
		tumor              // channel: [mandatory] [ val(group), val(meta), file(tumor_merged_vcf) ]
		normal             // channel: [mandatory] [ val(group), val(meta), file(normal_merged_vcf) ]
		meta               // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
	main:
		ch_versions = Channel.empty()

		GENE_INTERSECT ( tumor.mix(normal), params.gene_gtf )
		COYOTE_SEGMENTS_JSON ( GENE_INTERSECT.out.intersected )
		COYOTE_SEGMENTS ( tumor.mix(normal) )
		ch_versions = ch_versions.mix(COYOTE_SEGMENTS.out.versions)

		MERGE_SEGMENTS ( COYOTE_SEGMENTS.out.filtered.groupTuple().view() )

	emit:
		segments 	= 	MERGE_SEGMENTS.out.merged	// channel: [ val(group), file(cn-segments.panel.merged.bed) ]
		versions    =   ch_versions 				// channel: [ file(versions) ]
}