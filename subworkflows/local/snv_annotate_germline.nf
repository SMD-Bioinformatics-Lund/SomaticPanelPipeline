#!/usr/bin/env nextflow

include { ANNOTATE_VEP  as NORMAL_VEP           } from '../../modules/local/filters/main'

workflow SNV_ANNOTATE {
    take: 
        normal_germline // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        meta            // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()
        normal_no_index = normal_germline.map {
            group, meta, vcf, index ->
            tuple(group, meta, vcf)
        }
        // NON-OPTIONAL, NEEDED BY COYOTE
        NORMAL_VEP { normal_no_index }
        ch_versions = ch_versions.mix(NORMAL_VEP.out.versions)

    emit:
        vep_vcf             =   NORMAL_VEP.out.vcf_vep                           // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        versions            =   ch_versions                                      // channel: [ file(versions) ]
}