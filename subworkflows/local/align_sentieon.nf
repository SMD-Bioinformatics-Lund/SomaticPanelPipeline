#!/usr/bin/env nextflow

include { BWA_UMI                } from '../../modules/local/sentieon/main'
include { MARKDUP                } from '../../modules/local/sentieon/main'
include { BQSR                   } from '../../modules/local/sentieon/main'
include { DUMMY_DEDUP_METRICS    } from '../../modules/local/sentieon/main'




workflow ALIGN_SENTIEON {
    take: 
        fastq_input         // channel: [mandatory] [ val(meta), path([ reads ]) ]
        alt_bam_path        // channel: [optional]  [ val(meta), path([bam,bai]) ]
        meta                // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        BWA_UMI ( fastq_input )
        ch_versions = ch_versions.mix(BWA_UMI.out.versions)

        MARKDUP ( BWA_UMI.out.bam_umi_markdup )
        ch_versions = ch_versions.mix(MARKDUP.out.versions)

        DUMMY_DEDUP_METRICS( alt_bam_path )

        dedup_metrics = MARKDUP.out.dedup_metrics.mix(DUMMY_DEDUP_METRICS.out.dedup_metrics)
        
        BQSR ( BWA_UMI.out.bam_umi.mix(alt_bam_path) )

        bam_umi = BWA_UMI.out.bam_umi.mix(alt_bam_path)
        ch_versions = ch_versions.mix(BQSR.out.versions)

    emit:
        bam_dedup               =   MARKDUP.out.bam_dedup                   // channel: [ val(group), val(meta), file(bam), file(bai) ] 
        bam_umi                 =   bam_umi                                 // channel: [ val(group), val(meta), file(bam), file(bai) ]         
        bqsr_table              =   BQSR.out.bqsr_table                     // channel: [ val(group), val(meta), file(bqsr_table) ]
        dedup_metrics           =   dedup_metrics                           // channel: [ val(group), val(meta), file(dedup_metrics) ]
        versions                =   ch_versions                             // channel: [ file(versions) ]
        
}