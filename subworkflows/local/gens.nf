#!/usr/bin/env nextflow

include { GENS_YAML                   } from '../../modules/local/gens/main'

workflow GENS {
    take: 
        merged_gens_v4   // channel: [ val(group), val(meta), file(baf.bed.gz), file(cov.bed.gz) ] -- from MERGE_GENS.out.merged_gens_for_v4
        loh_cat          // channel: [ val(group), val(meta), file(*.loh_cat.tsv) ] -- from BIOMARKERS.out.loh_cat (always emitted, placeholder if !params.loh)
    
    main:   
        ch_versions = Channel.empty()

        // --- debug views ---
        merged_gens_v4.view { "MERGED_GENS_V4: $it" }

        // group per sample baf/cov into one per sample group
        ch_gens_v4_grouped  = merged_gens_v4
            .groupTuple()
            .map{ group, meta_list, baf_files, cov_files -> 
                tuple(group, meta_list, baf_files, cov_files)
            }
        
        ch_gens_v4_grouped.view { "GROUPED: $it" }
        loh_cat.view { "LOH_CAT: $it" }
        
        ch_gens_yaml_input = ch_gens_v4_grouped
            .join(loh_cat)
        
        ch_gens_yaml_input.view { "JOIN_INPUT: $it" }

        GENS_YAML ( ch_gens_yaml_input )
        ch_versions = ch_versions.mix(GENS_YAML.out.versions)

    emit:
        gens_yaml   =   GENS_YAML.out.gens_yaml     // channel: [ val(group), file(*.gens_somatic.yaml) ]
        versions    =   ch_versions                 // channel: [ file(versions) ]
}