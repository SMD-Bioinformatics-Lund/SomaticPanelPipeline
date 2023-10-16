#!/usr/bin/env nextflow

include { CREATE_SNVPON               } from '../../modules/local/filters/main'

workflow CREATE_SNV_PON {
    take: 
        concat_vcfs

    main:
        ch_versions = Channel.empty()

        CREATE_SNVPON( concat_vcfs.groupTuple(by:[1]))
        ch_versions = ch_versions.mix(CREATE_SNVPON.out.versions)

    emit:
        versions    =   ch_versions
        snv_pon     =   CREATE_SNVPON.out.SNV_PON
}