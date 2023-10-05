#!/usr/bin/env nextflow

include { CNVKITREF                     } from '../../modules/local/cnvkit/main'
include { CNVKITREF as CNVKITREF_EXONS  } from '../../modules/local/cnvkit/main'
include { CNVKITREF as CNVKITREF_BB     } from '../../modules/local/cnvkit/main'

workflow CNVKITREFS {
	take:
		sample		   		// val(id), file(cram), file(crai), file(bai)
        cnvkit_name    		// val(name)

	main:
		name_ch = channel.of(cnvkit_name)
		sample_bam = sample.map{ val-> tuple(val[0], val[1], val[2] ) }
        CNVKITREF(
			name_ch.combine(sample_bam).groupTuple(),
			"full",
			params.regions_bed
		)
		if (params.cnvkit_split) {
			CNVKITREF_EXONS(
				name_ch.combine(sample_bam).groupTuple(),
				"exons",
				params.regions_bed_exons
			)
			CNVKITREF_BB(
				name_ch.combine(sample_bam).groupTuple(),
				"backbone",
				params.regions_bed_backbone
			)
		}
		
	emit:
		reference = CNVKITREF.out.cnvkit_ref
        


}