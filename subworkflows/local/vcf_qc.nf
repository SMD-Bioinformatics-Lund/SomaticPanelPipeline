#!/usr/bin/env nextflow

include { CONTAMINATION            } from '../../modules/local/qc/main'
include { PEDDY                    } from '../../modules/local/peddy/main'

workflow VCF_QC {
    take:        
        vep_vcf                    // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        tumor_germline             // channel: [ val(group), file(vcf), file(tbi) ]
        normal_germline             // channel: [ val(group), val(meta) path(vcf), path(index) ]
        meta_ch                    // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        vep_vcf_filtered = vep_vcf.filter { group, meta, vcf ->
            meta.sex != false
        }

        ped_ch_tumor = meta_ch.filter { group, meta ->
            meta.type == "T"
        }.join(tumor_germline)

        ped_ch_normal = meta_ch.filter { group, meta ->
            meta.type == "N"
        }.join(normal_germline, by:[0,1])

        PEDDY { ped_ch_tumor.mix(ped_ch_normal) }

        CONTAMINATION { vep_vcf }
        ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
        
    emit:
        qcdone                  =   CONTAMINATION.out.contamination_cdm                  // channel: [ tuple val(group), file("dist.txt"), file("sampleid.png") ]
        results                 =   CONTAMINATION.out.contamination_result_files         // channel: [ val(group), val(meta), file(cdm) ]

}