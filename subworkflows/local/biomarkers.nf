#!/usr/bin/env nextflow

include { ASCAT_252                   } from '../../modules/local/ASCAT/main'
include { ASCAT_30                    } from '../../modules/local/ASCAT/main'
include { CNVKIT2ASCAT                } from '../../modules/local/ASCAT/main'

workflow BIOMARKERS {
    take: 
        baflogr

    main:
        CNVKIT2ASCAT ( baflogr )
        ASCAT_252 ( CNVKIT2ASCAT.out.ascat_input )
        ASCAT_30 ( CNVKIT2ASCAT.out.ascat_input )

    emit:
        ascat = CNVKIT2ASCAT.out.ascat_input
        //hrdscore_ascat252 = ASCAT_252.out.hrdscore
        //hrdscore_ascat30 = ASCAT_30.out.hrdscore

}