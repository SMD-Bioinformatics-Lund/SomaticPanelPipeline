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
            .filter { group, meta -> meta.type == 'T' || meta.type == 'tumor' }
            .join(tumor_germline.map { group, germline_meta, vcf, tbi -> tuple(group, vcf, tbi) })
            .map { group, meta, vcf, tbi -> tuple(group, meta, vcf, tbi) }

        ped_ch_normal = meta_filtered
            .filter { group, meta -> meta.type == 'N' || meta.type == 'normal' }
            .join(normal_germline, by:[0,1])

        PEDDY ( ped_ch_tumor.mix(ped_ch_normal) )
        PEDDY2CDM ( PEDDY.out.peddy_files )

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
