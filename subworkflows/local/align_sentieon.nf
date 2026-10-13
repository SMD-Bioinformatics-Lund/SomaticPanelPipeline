#!/usr/bin/env nextflow

include { UMI_CONSENSUS          } from '../../modules/local/sentieon/main'
include { BWA_UMI                } from '../../modules/local/sentieon/main'
include { BWA_SORT               } from '../../modules/local/sentieon/main'
include { LOCUS_COLLECTOR        } from '../../modules/local/sentieon/main'
include { MARKDUP                } from '../../modules/local/sentieon/main'
include { BQSR_UMI               } from '../../modules/local/sentieon/main'


workflow ALIGN_SENTIEON {
    take:
        fastq_input         // channel: [mandatory] [ val(group), val(meta), [ reads ] ]
        alt_bam_path        // channel: [optional]  [ val(group), val(meta), file(bam), file(bai) ]
        meta                // channel: [mandatory] [ val(group), val(meta) ]

    main:
        ch_versions = Channel.empty()

        UMI_CONSENSUS ( fastq_input )
        ch_versions = ch_versions.mix(UMI_CONSENSUS.out.versions)

        BWA_UMI ( UMI_CONSENSUS.out.consensus_fastq )
        ch_versions = ch_versions.mix(BWA_UMI.out.versions)

        BWA_SORT ( UMI_CONSENSUS.out.noumi_sam )
        ch_versions = ch_versions.mix(BWA_SORT.out.versions)

        LOCUS_COLLECTOR ( BWA_SORT.out.sorted_bam )
        ch_versions = ch_versions.mix(LOCUS_COLLECTOR.out.versions)

        MARKDUP ( LOCUS_COLLECTOR.out.score )
        ch_versions = ch_versions.mix(MARKDUP.out.versions)

        BQSR_UMI ( BWA_UMI.out.bam_umi.mix(alt_bam_path) )
        ch_versions = ch_versions.mix(BQSR_UMI.out.versions)

    emit:
        bam_dedup               =   MARKDUP.out.dedup_bam                   // channel: [ val(group), val(meta), file(bam), file(bai) ]
        bam_umi                 =   BQSR_UMI.out.bam_bqsr                   // channel: [ val(group), val(meta), file(bam), file(bai), file(bqsr.table) ]
        dedup_metrics           =   MARKDUP.out.metrics                     // channel: [ val(group), val(meta), file(dedup_metrics) ]
        versions                =   ch_versions                             // channel: [ file(versions) ]
        
}
