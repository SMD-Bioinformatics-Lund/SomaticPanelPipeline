#!/usr/bin/env nextflow

include { BWA_UMI                } from '../../modules/local/sentieon/main'
include { MARKDUP                } from '../../modules/local/sentieon/main'
include { BQSR_UMI               } from '../../modules/local/sentieon/main'


workflow ALIGN_SENTIEON {
    take: 
        fastq_input         // channel: [mandatory] [ val(meta), [ reads ] ]
        meta                // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        BWA_UMI ( fastq_input )
        ch_versions = ch_versions.mix(BWA_UMI.out.versions)

        MARKDUP ( BWA_UMI.out.bam_umi_markdup )
        ch_versions = ch_versions.mix(MARKDUP.out.versions)

        BQSR_UMI ( BWA_UMI.out.bam_umi )
        ch_versions = ch_versions.mix(BQSR_UMI.out.versions)

    emit:
        bam_dedup               =   MARKDUP.out.bam_qc                      // channel: [ val(group), val(meta), file(bam), file(bai) ]
        bam_umi                 =   BQSR_UMI.out.bam_varcall                // channel: [ val(group), val(meta), file(bam), file(bai), file(bqsr.table) ]   
        dedup_metrics           =   MARKDUP.out.dedup_metrics               // channel: [ val(group), val(meta), file(dedup_metrics) ]
        versions                =   ch_versions                             // channel: [ file(versions) ]
        
}