#!/usr/bin/env nextflow

include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { LOWCOV               } from '../../modules/local/qc/main'

workflow QC {
    take: 
        qc
        bam_lowcov

    main:

        QC_TO_CDM ( qc )
        LOWCOV ( bam_lowcov )

    emit:
        qcdone = QC_TO_CDM.out.cdm_done
        lowcov = LOWCOV.out.lowcov_regions

}