#!/usr/bin/env nextflow

include { CREATE_SNVPON               } from '../../modules/local/pon/main'

workflow CREATE_SNV_PON {
    take: 
        concat_vcfs     // channel: [group, meta, caller, vcf.gz, index]

    main:
        ch_versions = Channel.empty()

        CREATE_SNVPON( concat_vcfs.groupTuple(by:[2]))
        ch_versions = ch_versions.mix(CREATE_SNVPON.out.versions)

    emit:
        versions    =   ch_versions                     // channel: [ file(versions) ]
        snv_pon     =   CREATE_SNVPON.out.SNV_PON       // channel: [groups, metas, caller, PON.snv]
}
