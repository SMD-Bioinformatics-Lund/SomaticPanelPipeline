#!/usr/bin/env nextflow

include { GENS_YAML                   } from '../../modules/local/gens/main'

workflow GENS {
    take: 
        ouputs_for_new_gens_yaml   // channel: [ val(group), val(meta), file("*.final.sorted.baf.bed.gz*"), file("*.final.sorted.cov.bed.gz*")] -- from MERGE_GENS.out.ouputs_for_new_gens_yaml
        loh_cat         // channel: [ val(group), val(meta), val(part), file("*.loh_cat.tsv") ] -- from BIOMARKERS.out.loh_cat
    
    main:   
        ch_versions = Channel.empty()

        // --- debug views ---
        new_gens_yaml.view { "new_gens_yaml: $it" }

        // group per sample baf/cov into one per sample group
        ch_new_gens_yaml_grouped  = new_gens_yaml
            .groupTuple()
            .map{ group, meta_list, baf_files, cov_files -> 
                tuple(group, meta_list, baf_files, cov_files) // should be grouped by group and meta here?
            }
        
        ch_new_gens_yaml_grouped.view { "ch_new_gens_yaml_grouped: $it" }
        loh_cat.view { "LOH_CAT: $it" }
        

        loh_cat_ch = loh_cat.map { group, meta, part, tsv -> tuple(group, meta, tsv) }

        ch_gens_yaml_input = ch_new_gens_yaml_grouped
            .join(loh_cat_ch)
        
        ch_gens_yaml_input.view { "ch_gens_yaml_input: $it" }

        GENS_YAML ( ch_gens_yaml_input )

        ch_versions = ch_versions.mix(GENS_YAML.out.versions)

    emit:
        gens_yaml   =   GENS_YAML.out.gens_yaml     // channel: [ val(group), file(*.gens_somatic.yaml) ]
        versions    =   ch_versions                 // channel: [ file(versions) ]
}