#!/usr/bin/env nextflow

include { COYOTE_YAML          } from '../../modules/local/coyote/main'
include { COYOTE               } from '../../modules/local/coyote/main'
include { OUTPUT_FILES         } from '../../modules/local/coyote/main' // newly added process to ease data import into Coyote3
include { OUTPUTS_YAML_COYOTE  } from '../../modules/local/coyote/main' // newly added process to ease data import into Coyote3

workflow ADD_TO_DB {

    take: 
        vcf             // channel: [mandatory] [ val(group), val(meta), file(vcf) ]
        lowcov          // channel: [mandatory] [ val(group), val(meta.type), file(lowcov) ]
        lowcov_d4       // channel: [mandatory] [ val(group), val(meta.type), file(lowcov) ]
        segments        // channel: [optional] [ val(group), file(segments) ]
        s_json          // channel: [optional] [ val(group), file(segments) ]
        gens            // channel: [optional] [ val(group), val(meta), file(gens) ]
        gatcov_plot     // channel: [optional] [ val(group), file(plot) ]
        fusions         // channel: [optional] [ val(group), file(vcf) ]
        biomarkers      // channel: [optional] [ val(group), file(json) ]
        cnvkit_plot     // channel: [optional] [ val(group), val(meta), val(part), file(cnvkit_overview.png) ]

    main:
        ch_coyote_info = Channel.empty() // Gather data for json_INFO listing output files for export into Coyote3
        lc = lowcov.map{ val-> tuple(val[0], val[2] ) }
        lc_d4 = lowcov_d4.map{ val-> tuple(val[0], val[2] ) }
        cnvkit_plot.groupTuple().set { cnvkit_plot_ch }
        cnvkit_plot_ch.map { tuple ->
            if (tuple[1].size() >= 2) {
                int idx = tuple[1].findIndexOf { it['type'] == 'T' || it['type'] == 'tumor' }
                [tuple[0], tuple[3][idx]]
            } else {
                [tuple[0], tuple[3][0]]
            }
        }.set { cnvkit_png }

        ch_coyote_info = ch_coyote_info
        .mix(lowcov_d4
            .filter { group, type, file -> type == 'T' || type == 'tumor' }
            .map {group, type, file -> tuple(group, 'cov', file) }) // using directly the naming in coyote yaml for coyote_cov_json
        .mix(gatcov_plot
            .map { group, file -> tuple(group, 'cnvprofile_gatk', file) }) // file is modeled.png (preferred for coyote)
        .mix(cnvkit_plot
            .filter { group, meta, part, file -> meta.type == 'T' || meta.type == 'tumor'} 
            .map { group, meta, part, file -> tuple(group, 'cnvprofile_cnvkit', file) }) // file is cnvkit_overview.png (if modeled.png absent)
        .mix(s_json
            .map { group, file -> tuple(group, 'cnv', file) }) // using directly the naming in coyote yaml)
        .mix(fusions
            .map { group, file -> tuple(group, 'transloc', file) }) // using directly the naming in coyote yaml for merged.annotated.vcf
        .mix(biomarkers
            .map { group, file -> tuple(group, 'biomarkers', file) }) // using directly the naming in coyote yaml for bio.json
    
        optional = lc.mix(s_json,gatcov_plot,biomarkers,fusions,cnvkit_png).groupTuple()
        optional_json = lc_d4.mix(s_json,gatcov_plot,biomarkers,fusions,cnvkit_png).groupTuple()
        COYOTE { vcf.join(optional) }
        // COYOTE_YAML { vcf.join(optional_json) }
        ch_coyote_info_grouped = ch_coyote_info.groupTuple()
        OUTPUT_FILES ( ch_coyote_info_grouped ) // creates a json_INFO from the channel ch_coyote_info_grouped
        OUTPUTS_YAML_COYOTE ( vcf.join(OUTPUT_FILES.out.json_INFO) ) 

    emit:
        coyotedone = COYOTE.out.coyote_import        // channel: [ val(group), file(coyote) ]
        //coyotedone = COYOTE_YAML.out.coyote_import   // channel: [ val(group), file(coyote) ], old way to make the YAML for Coyote3
        coyoteyaml = OUTPUTS_YAML_COYOTE.out.yaml_coyote3 // channel: [ val(group), file(coyote) ]
}

