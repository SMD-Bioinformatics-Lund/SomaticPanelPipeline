#!/usr/bin/env nextflow

include { MSISENSOR_MSI               } from '../../modules/local/msisensor/main'
include { MSISENSOR_PRO               } from '../../modules/local/msisensor/main'
include { CNVKIT2SCARHRD              } from '../../modules/local/hrdsw/main'
include { SCARHRD                     } from '../../modules/local/hrdsw/main'
include { BIOMARKERS_TO_JSON          } from '../../modules/local/filters/main'


workflow BIOMARKERS {
    take:
        // meta                   // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        cnvkitsegments         // channel: [mandatory] [ val(group), val(meta), val(part(backbone)), file("${group}.${meta.id}.${part}.call*.cns") ]
        bam_umi                // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr) ]
        bam_dedup              // channel: [mandatory] [ val(group), val(meta), file(marked_bam), file(marked_bai), file(bqsr) ]

    main:
        ch_versions = Channel.empty()

        // HRD //
        CNVKIT2SCARHRD ( cnvkitsegments.filter { it -> it[1].type == "T" })
        ch_versions = ch_versions.mix(CNVKIT2SCARHRD.out.versions)

        SCARHRD ( CNVKIT2SCARHRD.out.scarHRD_segments)
        ch_versions = ch_versions.mix(SCARHRD.out.versions)

        // MSI //
        MSISENSOR_MSI (bam_umi.groupTuple())
        ch_versions = ch_versions.mix(MSISENSOR_MSI.out.versions)

        MSISENSOR_PRO (bam_umi.groupTuple())
        ch_versions = ch_versions.mix(MSISENSOR_PRO.out.versions)

        // combine biomarkers //
        BIOMARKERS_TO_JSON( SCARHRD.out.scarHRD_score.mix(MSISENSOR_PRO.out.msi_score,MSISENSOR_MSI.out.msi_score_paired).groupTuple() )
        ch_versions = ch_versions.mix(BIOMARKERS_TO_JSON.out.versions)

        output = Channel.empty().mix(BIOMARKERS_TO_JSON.out.biomarkers_json)

    emit:
        biomarkers  =   output          // channel: [ val(group), file(bio.json) ]
        versions    =   ch_versions     // channel: [ file(versions) ]

}
