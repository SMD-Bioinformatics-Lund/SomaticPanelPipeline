#!/usr/bin/env nextflow

include { MSISENSOR                   } from '../../modules/local/msisensor/main'
include { CNVKIT2OVAHRDSCAR           } from '../../modules/local/hrdsw/main'
include { CNVKIT2SCARHRD              } from '../../modules/local/hrdsw/main'
include { SCARHRD                     } from '../../modules/local/hrdsw/main'
include { BIOMARKERS_TO_JSON          } from '../../modules/local/filters/main'
include { CALL_LOH                     } from '../../modules/local/loh/main'


workflow BIOMARKERS {
    take: 
        meta                   // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        cnvkitsegments         // channel: [mandatory] [ val(group), val(meta), val(part), file("${group}.${meta.id}.${part}.call*.cns") ]
        cnvkitsegment_for_loh  // channel: [mandatory] [ val(group), val(meta), val(part), file("${group}.${meta.id}.${part}.call*.cns") ]
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


        // Loss of Heterozygosity (LOH) calculation and type evaluation
        // input for CALL_LOH is an output from the CNV_CALLING subworkflow and is cnvkitsegments from CNVKIT_CALL using "full" .cns (gmshem) 
        // or "backbone" .cns (solid).
        // params.loh in the process config in biomarkers.config
        cnvkitsegments_tumor = cnvkitsegment_for_loh.filter { group, meta, part, cns -> meta.type == "T" } 
        CALL_LOH ( cnvkitsegments_tumor )
        ch_versions = ch_versions.mix(CALL_LOH.out.versions)
 
        loh_bed_ch = Channel.empty().mix(CALL_LOH.out.loh_bed) // ensure channel is emitted even if params.loh is false and CALL_LOH not executed


    emit:
        biomarkers  =   output         // channel: [ val(group), file(bio.json) ]
        versions    =   ch_versions    // channel: [ file(versions) ]
        loh_bed     =   loh_bed_ch     // channel: [ val(group), val(meta, type="T"), val(part), file("*.loh_cat.bed") ]

}