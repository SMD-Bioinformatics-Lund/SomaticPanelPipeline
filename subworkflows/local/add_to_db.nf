#!/usr/bin/env nextflow

include { COYOTE_YAML          } from '../../modules/local/coyote/main'
include { COYOTE               } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf             // channel: [mandatory] [ val(group), val(meta), file(vcf) ]
        lowcov          // channel: [mandatory] [ val(group), val(meta), val(meta.type), file(lowcov) ]
        lowcov_d4       // channel: [mandatory] [ val(group), val(meta), val(meta.type), file(lowcov) ]
        segments        // channel: [optional] [ val(group), file(segments) ]
        s_json          // channel: [optional] [ val(group), file(segments) ]
        gens            // channel: [optional] [ val(group), val(meta), file(gens) ]
        gatcov_plot     // channel: [optional] [ val(group), file(plot) ]
        fusions         // channel: [optional] [ val(group), file(vcf) ]
        biomarkers      // channel: [optional] [ val(group), file(json) ]
        cnvkit_plot     // channel: [optional] [ val(group), val(meta), val(part), file(cnvkit_overview.png) ]

    main:
        lc = lowcov.map { row ->
            tuple(row[0], 'lowcov', row[3].name)
        }
        lc_d4 = lowcov_d4.map { row ->
            tuple(row[0], 'cov', row[3].name)
        }
        cnvkit_plot.groupTuple().set { cnvkit_plot_ch }
        cnvkit_plot_ch.map { row ->
            if (row[1].size() >= 2) {
                int idx = row[1].findIndexOf { it['type'] == 'T' || it['type'] == 'tumor' }
                tuple(row[0], 'cnvprofile', row[3][idx].name)
            } else {
                tuple(row[0], 'cnvprofile', row[3][0].name)
            }
        }.set { cnvkit_png }

        cnv_json = s_json.map { row ->
            tuple(row[0], 'cnv', row[-1].name)
        }
        gatk_plot = gatcov_plot.map { row ->
            tuple(row[0], 'cnvprofile', row[-1].name)
        }
        biomarkers_json = biomarkers.map { row ->
            tuple(row[0], 'biomarkers', row[-1].name)
        }
        fusion_vcf = fusions.map { row ->
            tuple(row[0], 'transloc', row[-1].name)
        }

        optional = lc.mix(cnv_json, gatk_plot, biomarkers_json, fusion_vcf, cnvkit_png).groupTuple(by:[0])
        optional_json = lc_d4.mix(cnv_json, gatk_plot, biomarkers_json, fusion_vcf, cnvkit_png).groupTuple(by:[0])

        coyote_input = vcf.join(optional).map { group, meta, vcf_file, import_types, import_files ->
            tuple(group, meta, vcf_file, import_types, import_files)
        }
        coyote_yaml_input = vcf.join(optional_json).map { group, meta, vcf_file, import_types, import_files ->
            tuple(group, meta, vcf_file, import_types, import_files)
        }

        COYOTE ( coyote_input )
        COYOTE_YAML ( coyote_yaml_input )
        coyote_done = COYOTE.out.coyote_import.mix(COYOTE_YAML.out.coyote_import)
        

    emit:
        coyotedone = coyote_done        // channel: [ val(group), file(coyote command or coyote yaml) ]
}
