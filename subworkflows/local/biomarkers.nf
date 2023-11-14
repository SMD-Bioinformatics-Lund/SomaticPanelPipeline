#!/usr/bin/env nextflow

include { MSISENSOR                   } from '../../modules/local/msisensor/main'
include { CNVKIT2OVAHRDSCAR           } from '../../modules/local/hrdsw/main'
include { CNVKIT2SCARHRD              } from '../../modules/local/hrdsw/main'
include { SCARHRD                     } from '../../modules/local/hrdsw/main'
include { BIOMARKERS_TO_JSON          } from '../../modules/local/filters/main'


workflow BIOMARKERS {
    take: 
        meta                   // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        cnvkitsegments         // channel: [mandatory] [ val(group), val(meta), val(part(backbone)), file("${group}.${meta.id}.${part}.call*.cns") ]
        bam_umi                // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr) ]
        bam_dedup              // channel: [mandatory] [ val(group), val(meta), file(marked_bam), file(marked_bai), file(bqsr) ]

    main:
        ch_versions = Channel.empty()

        if (params.other_biomarkers) {
            // HRD //
            if (params.hrd) {
                CNVKIT2SCARHRD ( cnvkitsegments.filter { it -> it[1].type == "T" })
                ch_versions = ch_versions.mix(CNVKIT2SCARHRD.out.versions)

                SCARHRD ( CNVKIT2SCARHRD.out.scarHRD_segments) 
                ch_versions = ch_versions.mix(SCARHRD.out.versions)
            }

            // MSI //
            if (params.msi) {
                MSISENSOR (bam_umi.groupTuple())
                ch_versions = ch_versions.mix(MSISENSOR.out.versions)
            }

            // combine biomarkers //
            BIOMARKERS_TO_JSON( SCARHRD.out.scarHRD_score.mix(MSISENSOR.out.msi_score,MSISENSOR.out.msi_score_paired).groupTuple() )
            ch_versions = ch_versions.mix(BIOMARKERS_TO_JSON.out.versions)

            output = BIOMARKERS_TO_JSON.out.biomarkers_json
        }
        else {
            output = Channel.empty()
        }

    emit:
        biomarkers  =   output          // channel: [ val(group), file(bio.json) ]
        versions    =   ch_versions     // channel: [ file(versions) ]

}