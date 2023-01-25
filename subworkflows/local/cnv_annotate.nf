#!/usr/bin/env nextflow

include { ANNOTATE_VEP                          } from '../../modules/local/filters/main'
include { COYOTE_SEGMENTS                       } from '../../modules/local/filters/main'
include { MERGE_SEGMENTS                        } from '../../modules/local/filters/main'


workflow CNV_ANNOTATE {
	take: 
		tumor              // val(group), val(meta), file(merged_vcf) // tumor
		normal             // val(group), val(meta), file(merged_vcf) // normal optional
		meta               // map: (csv meta info)

	main:
        ANNOTATE_VEP ( tumor.mix(normal) )
		// choose panel, combine with out from vep //
		COYOTE_SEGMENTS ( ANNOTATE_VEP.out.vcf_vep )
		MERGE_SEGMENTS ( COYOTE_SEGMENTS.out.filtered.groupTuple().view() )

	emit:
        annotated = ANNOTATE_VEP.out.vcf_vep
		segments = MERGE_SEGMENTS.out.merged

}