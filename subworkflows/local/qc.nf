#!/usr/bin/env nextflow

include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { LOWCOV               } from '../../modules/local/qc/main'
include { QC_VALUES            } from '../../modules/local/qc/main'

workflow QC {
    take: 
        qc              // channel: [ tuple val(group), val(meta), file("json.QC") ]
        bam_lowcov      // channel: [ val(group), val(meta), file(bam), file(bai), file(dedup_metrics.txt) ]

    main:
        ch_versions = Channel.empty()

        QC_VALUES ( qc )

        QC_TO_CDM ( qc )

        LOWCOV ( bam_lowcov )
        ch_versions = ch_versions.mix(LOWCOV.out.versions)

    emit:
        qcdone      =   QC_TO_CDM.out.cdm_done      // channel: [ val(group), val(meta), file(cdm) ]
        lowcov      =   LOWCOV.out.lowcov_regions   // channel: [ val(group), val(meta.type), file(lowcov.bed) ]
        melt_qc     =   QC_VALUES.out.qc_melt_val   // channel: [ val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) ]
        versions    =   ch_versions                 // channel: [ file(versions) ]
}