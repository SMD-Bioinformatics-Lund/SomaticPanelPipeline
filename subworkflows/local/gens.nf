#!/usr/bin/env nextflow

include { GENS_YAML                   } from '../../modules/local/gens/main'

workflow GENS {
    take: 
        ouputs_for_new_gens_yaml    // channel: [ val(group), val(meta), file("*.final.sorted.baf.bed.gz*"), file("*.final.sorted.cov.bed.gz*")] 
                                    // -- from MERGE_GENS.out.ouputs_for_new_gens_yaml
        loh_bed                     // channel: [ val(group), val(meta, type="T"), val(part), file("*.loh_cat.bed") ] 
                                    // -- one item per group (tumor only), from BIOMARKERS.out.loh_bed
    
    main:   
        ch_versions = Channel.empty()

        // Collapse per-sample rows into one per group:
        // → [ group, [meta_T, meta_N], [baf_T, baf_N], [cov_T, cov_N] ]
        ch_grouped  = ouputs_for_new_gens_yaml
            .groupTuple()
        
        if( params.loh ) {

            ch_gens_yaml_input = ch_grouped
                .join(loh_bed
                    .map { group, meta, part, bed ->
                        tuple(group, bed) // drop meta and part, not needed for join
                    }
                )
            } 
        else {

            ch_gens_yaml_input = ch_grouped
                .map { group, meta_list, baf_files, cov_files ->
                    tuple(group, meta_list, baf_files, cov_files, null)
                }
        }

        GENS_YAML ( ch_gens_yaml_input )

        ch_versions = ch_versions.mix(GENS_YAML.out.versions)

    emit:
        gens_yaml   =   GENS_YAML.out.gens_yaml     // channel: [ val(group), file(*.gens_somatic.yaml) ]
        versions    =   ch_versions                 // channel: [ file(versions) ]
}