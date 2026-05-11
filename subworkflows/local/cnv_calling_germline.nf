#!/usr/bin/env nextflow


include { GATK_COUNT_GERMLINE                  } from '../../modules/local/GATK/main'
include { GATK_CALL_PLOIDY                     } from '../../modules/local/GATK/main'
include { GATK_CALL_GERMLINE_CNV               } from '../../modules/local/GATK/main'
include { FILTER_MERGE_GATK                    } from '../../modules/local/GATK/main'
include { POSTPROCESS                          } from '../../modules/local/GATK/main'
include { SVDB_MERGE_PANEL as JOIN_NORMAL      } from '../../modules/local/svdb/main'
include { FILTER_MANTA as FILTER_MANTA_NORMAL  } from '../../modules/local/filters/main'


workflow CNV_CALLING_GERMLINE {
    take: 
        bam_umi              // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]
        meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        gatk_ref             // channel: [mandatory] [ val(interger), val(part_of_genome) used for germline gatk-calling ]

    main:
        ch_versions = Channel.empty()

        // Do germline calling for normal
        GATK_COUNT_GERMLINE ( bam_umi.filter { it -> it[1].type == "N" })
        ch_versions = ch_versions.mix(GATK_COUNT_GERMLINE.out.versions)

        GATK_CALL_PLOIDY ( GATK_COUNT_GERMLINE.out.count_germline )
        ch_versions = ch_versions.mix(GATK_CALL_PLOIDY.out.versions)

        GATK_CALL_GERMLINE_CNV( GATK_COUNT_GERMLINE.out.count_germline.join(GATK_CALL_PLOIDY.out.gatk_ploidy,by:[0,1]).groupTuple(by:[0,1]).combine(gatk_ref) )
        ch_versions = ch_versions.mix(GATK_CALL_GERMLINE_CNV.out.versions)

        CALLED = GATK_CALL_GERMLINE_CNV.out.gatk_call_germline.groupTuple(by:[0,1])
        PLOIDY = GATK_CALL_PLOIDY.out.gatk_ploidy.groupTuple(by:[0,1])

        POSTPROCESS ( CALLED.join(PLOIDY,by:[0,1]).combine(gatk_ref.groupTuple(by:[3])))
        ch_versions = ch_versions.mix(POSTPROCESS.out.versions)

        FILTER_MERGE_GATK ( POSTPROCESS.out.gatk_germline_segmentsvcf )
        ch_versions = ch_versions.mix(FILTER_MERGE_GATK.out.versions)

        GATK_NORMAL = FILTER_MERGE_GATK.out.gatk_normal_vcf.join(meta.filter( it -> it[1].type == "N" ) ).map{ val-> tuple(val[0], val[2], val[1] )}
        JOIN_NORMAL ( GATK_NORMAL.groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_NORMAL.out.versions)


    emit:
        normal_vcf  =   JOIN_NORMAL.out.merged_vcf      // channel: [ val(group), val(vc), file(normal.merged.vcf) ]
        versions    =   ch_versions                     // channel: [ file(versions) ]

}