#!/usr/bin/env nextflow

include { VALIDATE_COYOTE_SNV            } from '../../modules/local/validate/main'

workflow SNV_VALIDATE {
    take:        
        final_vcf                    // channel: [ val(group), val(meta), file("final.vcf") ]
        assay_config                 // channel : [ file(params.assay_config) ]
        known_variants               // channel: [ file(params.known_validation_snvs) ]

    main:
        ch_versions = Channel.empty()

        VALIDATE_COYOTE_SNV ( final_vcf,assay_config,known_variants )
        
    emit:
	placeholder = VALIDATE_COYOTE_SNV.out.results

}
