#!/usr/bin/env nextflow

include { ANNOTATE_VEP                          } from '../../modules/local/filters/main'


workflow CNV_ANNOTATE {
	take: 
		tumor              // val(group), val(meta), file(merged_vcf) // tumor
		normal             // val(group), val(meta), file(merged_vcf) // normal optional
		meta               // map: (csv meta info)

	main:
        ANNOTATE_VEP ( tumor.mix(normal) )
		

	emit:
        annotated = ANNOTATE_VEP.out.vcf_vep

}