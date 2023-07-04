#!/usr/bin/env nextflow

include { MSISENSOR                   } from '../../modules/local/msisensor/main'
include { CNVKIT2OVAHRDSCAR           } from '../../modules/local/hrdsw/main'
include { CNVKIT2SCARHRD              } from '../../modules/local/hrdsw/main'
include { SCARHRD                     } from '../../modules/local/hrdsw/main'
include { BIOMARKERS_TO_JSON          } from '../../modules/local/filters/main'


workflow BIOMARKERS {
    take: 
        meta                   // meta info for sample(s)
        cnvkitsegments         // val(group), val(meta), val(part(backbone)), file("${group}.${meta.id}.${part}.call*.cns")
        bam_umi                // val(group), val(meta), file(bam), file(bai), file(bqsr) umi deduped bam
        bam_dedup              // val(group), val(meta), file(bam), file(bai), file(bqsr) markdup bam

    main:
        if (params.other_biomarkers) {
            // HRD //
            if (params.hrd) {
                CNVKIT2SCARHRD ( cnvkitsegments.filter { it -> it[1].type == "T" })
                SCARHRD ( CNVKIT2SCARHRD.out.scarHRD_segments) 
            }
            // MSI //
            if (params.msi) {
                MSISENSOR(bam_umi.groupTuple())
            }

            // combine biomarkers //
            BIOMARKERS_TO_JSON( SCARHRD.out.scarHRD_score.mix(MSISENSOR.out.msi_score,MSISENSOR.out.msi_score_paired).groupTuple() )

            output = BIOMARKERS_TO_JSON.out.biomarkers_json
        }
        else {
            output = Channel.empty()
        }




        // OLD ALTERNATE HRD ALGORITHMS AND SEGMENTATIONS //
        // CNVKIT2ASCAT ( baflogr )
        // //ASCAT_252 ( CNVKIT2ASCAT.out.ascat_input )
        // ASCAT_30 ( CNVKIT2ASCAT.out.ascat_input )
        // CNVKIT2OVAHRDSCAR ( cnvkitsegments.mix(cnvkitsegment_purity) )
        // ASCAT2SCARHRD ( ASCAT_30.out.baflogr ) //.join(ASCAT_30.out.ploidy, by:[0,1])
        // ASCAT2OVAHRDSCAR ( ASCAT_30.out.baflogr )
        // OVAHRDSCAR ( ASCAT2OVAHRDSCAR.out.ovaHRDscar_segments.mix(CNVKIT2OVAHRDSCAR.out.ovaHRDscar_segments) )

    emit:
        biomarkers = output

}