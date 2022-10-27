#!/usr/bin/env nextflow

include { QC_TO_CDM            } from '../../modules/local/qc/main'


workflow QC {
    take: 
        qc
        meta_QC

    main:

        QC_TO_CDM ( qc.join(meta_QC) )

    emit:
        qcdone = QC_TO_CDM.out.cdm_done



}