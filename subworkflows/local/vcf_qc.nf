#!/usr/bin/env nextflow

include { CONTAMINATION            } from '../../modules/local/qc/main'
include { CONTAMINATION_CDM        } from '../../modules/local/qc/main'
include { PEDDY                    } from '../../modules/local/peddy/main'
include { PEDDY2CDM                } from '../../modules/local/qc/main'

workflow VCF_QC {
    take:        
        vep_vcf                    // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        tumor_germline             // channel: [ val(group), val(meta), file(vcf), file(tbi) ]
        normal_germline            // channel: [ val(group), val(meta), file(vcf), file(tbi) ]
        meta_ch                    // channel: [ val(group), val(meta) ]

    main:
        ch_versions = Channel.empty()

        meta_filtered = meta_ch.filter { group, meta ->
            meta.sex != false
        }

        ped_ch_tumor = meta_filtered
            .filter { group, meta -> SampleGroups.isTumor(meta) }
            .map { group, meta -> tuple(group, meta.id, meta) }
            .join(tumor_germline.map { group, germline_meta, vcf, tbi ->
                def tumor_meta = SampleGroups.validateMetadata(group, germline_meta).find { SampleGroups.isTumor(it) }
                tumor_meta ? tuple(group, tumor_meta.id, vcf, tbi) : null
            }, by: [0,1], failOnDuplicate: true)
            .map { group, id, meta, vcf, tbi -> tuple(group, meta, vcf, tbi) }

        ped_ch_normal = meta_filtered
            .filter { group, meta -> SampleGroups.isNormal(meta) }
            .map { group, meta -> tuple(group, meta.id, meta) }
            .join(normal_germline.map { group, meta, vcf, tbi -> tuple(group, meta.id, vcf, tbi) }, by: [0,1], failOnDuplicate: true)
            .map { group, id, meta, vcf, tbi -> tuple(group, meta, vcf, tbi) }

        PEDDY ( ped_ch_tumor.mix(ped_ch_normal) )
        PEDDY2CDM ( PEDDY.out.peddy_files )
        ch_versions = ch_versions.mix(PEDDY2CDM.out.versions)

        CONTAMINATION { vep_vcf }
        ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
        ch_versions = ch_versions.mix(PEDDY.out.versions)

        CONTAMINATION_CDM { CONTAMINATION.out.contamination_values }
        
    emit:
        qcdone                  =   CONTAMINATION_CDM.out.contamination_cdm.mix(PEDDY2CDM.out.cdm)
        results                 =   CONTAMINATION.out.contamination_result_files
        peddy                   =   PEDDY2CDM.out.json
        versions                =   ch_versions

}
