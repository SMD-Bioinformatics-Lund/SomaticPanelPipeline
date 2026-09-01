#!/usr/bin/env nextflow

include { CNVKIT                                                   } from './cnvkit_cnv_calling'
include { GATK                                                     } from './gatk_cnv_calling'
include { MANTA                                                     } from '../../modules/local/manta/main'
include { SVDB_MERGE_PANEL as JOIN_TUMOR                            } from '../../modules/local/svdb/main'
include { SVDB_MERGE_PANEL as JOIN_NORMAL                           } from '../../modules/local/svdb/main'
include { FILTER_MANTA as FILTER_MANTA_TUMOR                        } from '../../modules/local/filters/main'
include { FILTER_MANTA as FILTER_MANTA_NORMAL                       } from '../../modules/local/filters/main'


workflow CNV_CALLING {
    take: 
        bam_umi              // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]
        germline_variants    // channel: [mandatory] [ val(group), val(meta), file(vcf), file(tbi) ]
        bam_markdup          // channel: [mandatory] [ val(group), val(meta), file(dedup_bam), file(dedup_bai)]  ]
        gatk_ref             // channel: [mandatory] [ val(interger), val(part_of_genome) used for germline gatk-calling ]

    main:
        ch_versions = Channel.empty()

        CNVKIT ( bam_umi, germline_variants )
        ch_versions = ch_versions.mix(CNVKIT.out.versions)
        CNVKIT_VCF_TUMOR = CNVKIT.out.tumor_vcf

        GATK ( bam_umi, gatk_ref )
        ch_versions = ch_versions.mix(GATK.out.versions)

        /////////////////////////// MANTA /////////////////////////////////////////////////////
        MANTA ( bam_markdup.groupTuple(), params.bedgz, "CNV" )
        ch_versions = ch_versions.mix(MANTA.out.versions)

        MANTA_NORMAL = MANTA.out.manta_vcf.flatMap { group, meta, vcfs ->
            files = vcfs instanceof List ? vcfs : [vcfs]
            normal_meta = meta.find { it.type == "N" || it.type == "normal" }
            normal_meta ? files.findAll { it.name.startsWith("${normal_meta.id}_manta.") }.collect { vcf -> tuple(group, normal_meta, vcf) } : []
        }

        FILTER_MANTA_NORMAL(MANTA_NORMAL)
        ch_versions = ch_versions.mix(FILTER_MANTA_NORMAL.out.versions)

        GATK_NORMAL_LABELLED = GATK.out.normal_vcf.map { group, meta, vcf ->
            tuple(group, meta, 'gatk', vcf)
        }
        MANTA_NORMAL_LABELLED = FILTER_MANTA_NORMAL.out.filtered.map { group, meta, vcf ->
            tuple(group, meta, 'manta', vcf)
        }

        JOIN_NORMAL ( GATK_NORMAL_LABELLED.mix(MANTA_NORMAL_LABELLED).groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_NORMAL.out.versions)

        // Join tumor vcf
        GATK_TUMOR = GATK.out.tumor_vcf

        MANTA_TUMOR = MANTA.out.manta_vcf.flatMap { group, meta, vcfs ->
            files = vcfs instanceof List ? vcfs : [vcfs]
            tumor_meta = meta.find { it.type == "T" || it.type == "tumor" }
            tumor_meta ? files.findAll { it.name.startsWith("${tumor_meta.id}_manta.") }.collect { vcf -> tuple(group, tumor_meta, vcf) } : []
        }

        FILTER_MANTA_TUMOR(MANTA_TUMOR)
        ch_versions = ch_versions.mix(FILTER_MANTA_TUMOR.out.versions)

        GATK_TUMOR_LABELLED = GATK_TUMOR.map { group, meta, vcf ->
            tuple(group, meta, 'gatk', vcf)
        }
        MANTA_TUMOR_LABELLED = FILTER_MANTA_TUMOR.out.filtered.map { group, meta, vcf ->
            tuple(group, meta, 'manta', vcf)
        }
        CNVKIT_TUMOR_LABELLED = CNVKIT_VCF_TUMOR.map { group, meta, vcf ->
            tuple(group, meta, 'cnvkit', vcf)
        }

        JOIN_TUMOR ( GATK_TUMOR_LABELLED.mix(MANTA_TUMOR_LABELLED, CNVKIT_TUMOR_LABELLED).groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_TUMOR.out.versions)

    emit:
        gatcov_plot =   GATK.out.gatcov_plot            // channel: [ val(group), val(meta), file(modeled.png) ]
        cnvkit_plot =   CNVKIT.out.plot                 // channel: [ val(group), val(meta), val(part), file(cnvkit_overview.png) ]
        cnvkit_hrd  =   CNVKIT.out.hrd                  // channel: [ val(group), val(meta), val(part), file(call.cns) ]
        tumor_vcf   =   JOIN_TUMOR.out.merged_vcf       // channel: [ val(group), val(meta), file(tumor.merged.vcf) ]
        normal_vcf  =   JOIN_NORMAL.out.merged_vcf      // channel: [ val(group), val(meta), file(normal.merged.vcf) ]
        gens        =   CNVKIT.out.gens                 // channel: [ val(group), val(meta), file(gens) ]
        gens_v4     =   CNVKIT.out.gens_v4              // channel: [ val(group), val(meta), file(gens_v4) ]
        versions    =   ch_versions                     // channel: [ file(versions) ]

}
