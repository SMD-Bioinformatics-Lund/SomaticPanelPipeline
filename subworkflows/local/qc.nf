#!/usr/bin/env nextflow

include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { LOWCOV               } from '../../modules/local/qc/main'
include { QC_VALUES            } from '../../modules/local/qc/main'

workflow QC {
    take: 
        qc                  // tuple val(group), val(meta), file("json.QC")
        bam_lowcov  

    main:
        ch_versions = Channel.empty()

        QC_VALUES( qc )

        QC_TO_CDM ( qc )

        LOWCOV ( bam_lowcov )
        ch_versions = ch_versions.mix(LOWCOV.out.versions)

    emit:
        qcdone      =   QC_TO_CDM.out.cdm_done
        lowcov      =   LOWCOV.out.lowcov_regions
        melt_qc     =   QC_VALUES.out.qc_melt_val
        versions    =   ch_versions 
}