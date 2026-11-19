#!/usr/bin/env nextflow

include { ANNOTATE_VEP as ANNOTATE_VEP_GERMLINE } from '../../modules/local/vep/main'

workflow SNV_ANNOTATE_GERMLINE {
    take:
        normal_germline // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        meta            // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()
        normal_no_index = normal_germline.map {
            group, meta_data, vcf, index ->
            tuple(group, meta_data, vcf)
        }
        // NON-OPTIONAL, NEEDED BY COYOTE
        ANNOTATE_VEP_GERMLINE { normal_no_index }
        ch_versions = ch_versions.mix(ANNOTATE_VEP_GERMLINE.out.versions)

    emit:
        vep_vcf             =   ANNOTATE_VEP_GERMLINE.out.vcf_vep                           // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        versions            =   ch_versions                                      // channel: [ file(versions) ]
}
