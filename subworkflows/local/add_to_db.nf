#!/usr/bin/env nextflow

include { COYOTE              } from '../../modules/local/coyote/main'
include { OUTPUT_FILES        } from '../../modules/local/coyote/main'
include { OUTPUTS_YAML_COYOTE } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {

    take:
        vcf             // channel: [mandatory] [ val(group), val(meta), file(vcf) ]
        lowcov          // channel: [mandatory] [ val(group), val(meta), val(meta.type), file(lowcov) ]
        lowcov_d4       // channel: [mandatory] [ val(group), val(meta), val(meta.type), file(lowcov) ]
        segments        // channel: [optional] [ val(group), file(segments) ]
        s_json          // channel: [optional] [ val(group), file(segments) ]
        gens            // channel: [optional] [ val(group), val(meta), file(gens) ]
        gatcov_plot     // channel: [optional] [ val(group), val(meta), file(plot) ]
        fusions         // channel: [optional] [ val(group), file(vcf) ]
        biomarkers      // channel: [optional] [ val(group), file(json) ]
        cnvkit_plot     // channel: [optional] [ val(group), val(meta), val(part), file(cnvkit_overview.png) ]
        bam_umi         // channel: [optional] [ val(group), val(meta), file(bam), file(bai), file(bqsr) ]

    main:
        ch_coyote_vcf = vcf.filter { group, meta, vcf_file ->
            SampleGroups.validateMetadata(group, meta)
        }
        lowcov_for_coyote = lowcov.map { group, meta, type, lowcov_file ->
            tuple(group, 'lowcov', lowcov_file.name)
        }

        cnvkit_png_for_coyote = cnvkit_plot
            .filter { group, meta, part, plot -> SampleGroups.isTumor(meta) }
            .map { group, meta, part, plot -> tuple(group, 'cnvprofile_cnvkit', plot.name) }

        cnv_json_for_coyote = s_json.map { group, cnv_file ->
            tuple(group, 'cnv', cnv_file.name)
        }
        gatk_plot_for_coyote = gatcov_plot.map { group, meta, plot_file ->
            tuple(group, 'cnvprofile_gatk', plot_file.name)
        }
        biomarkers_for_coyote = biomarkers.map { group, biomarker_file ->
            tuple(group, 'biomarkers', biomarker_file.name)
        }
        fusions_for_coyote = fusions.map { group, fusion_file ->
            tuple(group, 'transloc', fusion_file.name)
        }

        coyote_optional = lowcov_for_coyote
            .mix(cnv_json_for_coyote, gatk_plot_for_coyote, biomarkers_for_coyote, fusions_for_coyote, cnvkit_png_for_coyote)
            .groupTuple(by: [0])

        coyote_input = ch_coyote_vcf.join(coyote_optional, remainder: true, failOnDuplicate: true)
            .filter { group, meta, vcf_file, import_types, import_files -> vcf_file != null }
            .map { group, meta, vcf_file, import_types, import_files ->
                tuple(group, meta, vcf_file, import_types ?: [], import_files ?: [])
            }

        coyote_cov_json = lowcov_d4
            .filter { group, meta, type, cov_file ->
                type.toString().toLowerCase() in ['t', 'tumor']
            }
            .map { group, meta, type, cov_file ->
                tuple(group, 'cov', cov_file)
            }
        coyote_gatk_plot = gatcov_plot.map { group, meta, plot_file ->
            tuple(group, 'cnvprofile_gatk', plot_file)
        }
        coyote_cnvkit_plot = cnvkit_plot
            .filter { group, meta, part, plot_file ->
                meta.type.toString().toLowerCase() in ['t', 'tumor']
            }
            .map { group, meta, part, plot_file ->
                tuple(group, 'cnvprofile_cnvkit', plot_file)
            }
        coyote_cnv_json = s_json.map { group, cnv_file ->
            tuple(group, 'cnv', cnv_file)
        }
        coyote_fusions = fusions.map { group, fusion_file ->
            tuple(group, 'transloc', fusion_file)
        }
        coyote_biomarkers = biomarkers.map { group, biomarker_file ->
            tuple(group, 'biomarkers', biomarker_file)
        }

        coyote_info = coyote_cov_json
            .mix(coyote_gatk_plot, coyote_cnvkit_plot, coyote_cnv_json, coyote_fusions, coyote_biomarkers)
            .groupTuple(by: [0])

        COYOTE ( coyote_input )
        coyote_bams = bam_umi
            .map { group, meta, bam, bai, bqsr ->
                if (!SampleGroups.validateMetadata(group, meta)) return null
                tuple(group, meta.id, [bam: bam.name, bai: bai.name])
            }
            .groupTuple()
            .map { group, ids, bam_files ->
                if (ids.unique(false).size() != ids.size()) {
                    SampleGroups.skip(group, 'duplicate BAM/BAI entries')
                    return null
                }
                tuple(group, [ids, bam_files].transpose().collectEntries { id, files -> [(id): files] })
            }

        coyote_output_files = ch_coyote_vcf.join(coyote_info, remainder: true, failOnDuplicate: true)
            .filter { group, meta, vcf_file, labels, files -> vcf_file != null }
            .join(coyote_bams, remainder: true, failOnDuplicate: true)
            .filter { group, meta, vcf_file, labels, files, bam_files -> vcf_file != null }
            .map { group, meta, vcf_file, labels, files, bam_files ->
                def samples = SampleGroups.validateMetadata(group, meta)
                if (!samples) return null
                def case_meta = samples.find { SampleGroups.isTumor(it) } ?: samples[0]
                def control_meta = samples.find { it.id != case_meta.id && SampleGroups.isNormal(it) }
                def case_files = bam_files?.get(case_meta.id)
                def control_files = control_meta ? bam_files?.get(control_meta.id) : null
                if ((case_files?.bam && case_files.bam == control_files?.bam) ||
                    (case_files?.bai && case_files.bai == control_files?.bai)) {
                    SampleGroups.skip(group, 'case and control BAM/BAI publication names collide')
                    return null
                }
                def sample_ids = samples.collect { it.id } as Set
                if (bam_files && !sample_ids.containsAll(bam_files.keySet())) {
                    SampleGroups.skip(group, 'BAM metadata contains foreign sample IDs')
                    return null
                }
                tuple(group, labels ?: [], files ?: [], [
                    case_bam: case_files?.bam,
                    case_bai: case_files?.bai,
                    control_bam: control_files?.bam,
                    control_bai: control_files?.bai
                ])
            }

        OUTPUT_FILES ( coyote_output_files )
        OUTPUTS_YAML_COYOTE ( ch_coyote_vcf.join(OUTPUT_FILES.out.json_INFO, by: 0, failOnDuplicate: true) )

    emit:
        coyotedone = COYOTE.out.coyote_import
        coyoteyaml = OUTPUTS_YAML_COYOTE.out.yaml_coyote3
}
