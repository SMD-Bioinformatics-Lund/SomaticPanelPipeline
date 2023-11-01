#!/usr/bin/env nextflow

include { CREATE_SNVPON               } from '../../modules/local/filters/main'

workflow CREATE_SNV_PON {
    take: 
        concat_vcfs     // channel: [mandatory] [ val(group), val(vc), file(vcf.gz) ]

    main:
        ch_versions = Channel.empty()

        CREATE_SNVPON( concat_vcfs.groupTuple(by:[1]))
        ch_versions = ch_versions.mix(CREATE_SNVPON.out.versions)

    emit:
        versions    =   ch_versions                     // channel: [ file(versions) ]
        snv_pon     =   CREATE_SNVPON.out.SNV_PON       // channel: [ val(group), val(vc), file(PON.snv) ]
}