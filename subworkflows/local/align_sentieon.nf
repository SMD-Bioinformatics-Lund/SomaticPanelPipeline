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
        ch_versions = Channel.empty()

        BWA_UMI ( fastq_input )
        ch_versions = ch_versions.mix(BWA_UMI.out.versions)

        MARKDUP ( BWA_UMI.out.bam_umi_markdup )
        ch_versions = ch_versions.mix(MARKDUP.out.versions)

        BQSR_UMI ( BWA_UMI.out.bam_umi )
        ch_versions = ch_versions.mix(BQSR_UMI.out.versions)

        SENTIEON_QC ( MARKDUP.out.bam_qc )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)

    emit:
        bam_lowcov  =   MARKDUP.out.bam_qc
        bam_umi     =   BQSR_UMI.out.bam_varcall
        qc_out      =   SENTIEON_QC.out.qc_cdm
        bam_dedup   =   MARKDUP.out.bam_bqsr
        versions    =   ch_versions 
}