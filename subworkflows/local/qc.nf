#!/usr/bin/env nextflow

include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { LOWCOV               } from '../../modules/local/qc/main'

workflow QC {
    take: 
        qc
        meta_QC
        bam_lowcov

    main:

        QC_TO_CDM ( qc.join(meta_QC) )
        LOWCOV ( bam_lowcov )

    emit:
        qcdone = QC_TO_CDM.out.cdm_done
        lowcov = LOWCOV.out.lowcov_regions

}