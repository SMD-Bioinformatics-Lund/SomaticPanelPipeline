#!/usr/bin/env nextflow

include { BWA_UMI                } from '../../modules/local/sentieon/main'
include { MARKDUP                } from '../../modules/local/sentieon/main'
include { BQSR_UMI               } from '../../modules/local/sentieon/main'



workflow ALIGN_SENTIEON {
    take: 
        fastq_input

    main:
        BWA_UMI ( fastq_input )
        MARKDUP (BWA_UMI.out.bam_umi)


    emit:
        bam = MARKDUP.out.bam_bqsr



}