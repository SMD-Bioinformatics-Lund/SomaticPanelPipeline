#!/usr/bin/env nextflow

include { GENS_YAML                   } from '../../modules/local/gens/main'

workflow GENS {
    take: 
        ouputs_for_new_gens_yaml    // channel: [ val(group), val(meta), file("*.final.sorted.baf.bed.gz*"), file("*.final.sorted.cov.bed.gz*")] 
                                    // -- from MERGE_GENS.out.ouputs_for_new_gens_yaml
        loh_cat                     // channel: [ val(group), val(meta, type="T"), val(part), file("*.loh_cat.tsv") ] 
                                    // -- one item per group (tumor only), from BIOMARKERS.out.loh_cat
    
    main:   
        ch_versions = Channel.empty()

        // Collapse per-sample rows into one per group:
        // → [ group, [meta_T, meta_N], [baf_T, baf_N], [cov_T, cov_N] ]
        ch_grouped  = ouputs_for_new_gens_yaml
            .groupTuple()
        
        // Drop meta and part for loh_cat - onyl group is needed to join with ch_grouped
        ch_loh = loh_cat
            .map { group, meta, part, tsv 
                -> tuple(group, tsv) }
        
        // Join on group → [ group, [meta_T, meta_N], [baf_T, baf_N], [cov_T, cov_N], tsv ]
        ch_gens_yaml_input = ch_grouped
            .join(ch_loh)

        GENS_YAML ( ch_gens_yaml_input )

        ch_versions = ch_versions.mix(GENS_YAML.out.versions)

    emit:
        gens_yaml   =   GENS_YAML.out.gens_yaml     // channel: [ val(group), file(*.gens_somatic.yaml) ]
        versions    =   ch_versions                 // channel: [ file(versions) ]
}