#!/usr/bin/env nextflow

include { BWA_UMI                } from '../../modules/local/sentieon/main'
include { MARKDUP                } from '../../modules/local/sentieon/main'
include { BQSR_UMI               } from '../../modules/local/sentieon/main'
include { SENTIEON_QC            } from '../../modules/local/sentieon/main'


workflow ALIGN_SENTIEON {
    take: 
        fastq_input
        meta

    main:
        BWA_UMI ( fastq_input )
        MARKDUP ( BWA_UMI.out.bam_umi_markdup )
        BQSR_UMI ( BWA_UMI.out.bam_umi )
        SENTIEON_QC ( MARKDUP.out.bam_qc )

    emit:
        bam_lowcov = MARKDUP.out.bam_qc
        bam_umi = BQSR_UMI.out.bam_varcall
        qc_out = SENTIEON_QC.out.qc_cdm
        bam_dedup = MARKDUP.out.bam_bqsr



}