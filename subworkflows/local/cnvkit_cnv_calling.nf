#!/usr/bin/env nextflow

include { CNVKIT_BATCH as CNVKIT_FULL                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_BACKBONE                           } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_EXONS                              } from '../../modules/local/cnvkit/main'
include { CNVKIT_GENS                                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_PLOT                                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL                                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_EXPORT_VCF                                         } from '../../modules/local/cnvkit/main'
include { CNVKIT_EXPORT_NEXUS_OGT                                   } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL as CNVKIT_CALL_TC                             } from '../../modules/local/cnvkit/main'
include { CNVKIT_EXPORT_VCF as CNVKIT_EXPORT_VCF_TC                 } from '../../modules/local/cnvkit/main'
include { MERGE_GENS                                                } from '../../modules/local/cnvkit/main'
include { TABIX_BGZIPTABIX as CNVKIT_BAF_BGZIPTABIX                 } from '../../modules/local/tabix/main'
include { TABIX_BGZIPTABIX as CNVKIT_COV_BGZIPTABIX                 } from '../../modules/local/tabix/main'


workflow CNVKIT {
    take:
        bam_umi              // channel: [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]
        germline_variants    // channel: [ val(group), val(meta), file(vcf), file(tbi) ]

    main:
        ch_versions = Channel.empty()

        if (!params.cnvkit_split) {
            CNVKIT_FULL ( bam_umi, params.cnvkit_reference, "full" )
            batch_plot_cns = CNVKIT_FULL.out.cnvkit_cns
            batch_plot_cnr = CNVKIT_FULL.out.cnvkit_cnr
            ch_versions = ch_versions.mix(CNVKIT_FULL.out.versions)

            cnvkit_full_input = batch_plot_cns.join(batch_plot_cnr, by:[0,1,3])
                .combine(germline_variants, by:[0])
                .map { group, meta, part, cns, cnr, germline_meta, vcf, tbi ->
                    tuple(group, meta, part, cns, cnr, vcf, tbi)
                }

            cnvkit_full_call_input = cnvkit_full_input.map { group, meta, part, cns, cnr, vcf, tbi ->
                tuple(group, meta, part, cns, vcf, tbi)
            }

            cnvkit_nexus_input = cnvkit_full_input.map { group, meta, part, cns, cnr, vcf, tbi ->
                tuple(group, meta, part, cnr, vcf, tbi)
            }

            CNVKIT_PLOT ( cnvkit_full_input )
            ch_versions = ch_versions.mix(CNVKIT_PLOT.out.versions)

            CNVKIT_CALL ( cnvkit_full_call_input, "false" )
            ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

            CNVKIT_EXPORT_VCF ( CNVKIT_CALL.out.cnvkitsegment )
            ch_versions = ch_versions.mix(CNVKIT_EXPORT_VCF.out.versions)

            CNVKIT_EXPORT_NEXUS_OGT ( cnvkit_nexus_input )
            ch_versions = ch_versions.mix(CNVKIT_EXPORT_NEXUS_OGT.out.versions)

            cnvkit_full_gens_input = batch_plot_cnr.combine(germline_variants, by:[0])
                .map { group, meta, cnr, part, germline_meta, vcf, tbi ->
                    tuple(group, meta, cnr, part, vcf, tbi)
                }

            CNVKIT_GENS ( cnvkit_full_gens_input )
            ch_versions = ch_versions.mix(CNVKIT_GENS.out.versions)

            ch_cnvkit_baf = CNVKIT_GENS.out.cnvkit_gens.map { group, meta, part, baf, cov ->
                tuple(group, meta, part, baf)
            }
            ch_cnvkit_cov = CNVKIT_GENS.out.cnvkit_gens.map { group, meta, part, baf, cov ->
                tuple(group, meta, part, cov)
            }
            CNVKIT_BAF_BGZIPTABIX ( ch_cnvkit_baf )
            CNVKIT_COV_BGZIPTABIX ( ch_cnvkit_cov )
            cnvkit_gens_compressed = CNVKIT_BAF_BGZIPTABIX.out.gz_index
                .join(CNVKIT_COV_BGZIPTABIX.out.gz_index, by:[0,1,2])
                .map { group, meta, part, baf_gz, baf_tbi, cov_gz, cov_tbi ->
                    tuple(group, meta, baf_gz, cov_gz)
                }

            MERGE_GENS ( cnvkit_gens_compressed )
            ch_versions = ch_versions.mix(CNVKIT_BAF_BGZIPTABIX.out.versions)
            ch_versions = ch_versions.mix(CNVKIT_COV_BGZIPTABIX.out.versions)
            ch_versions = ch_versions.mix(MERGE_GENS.out.versions)

            cnvkitplot = CNVKIT_PLOT.out.cnvkitplot
            cnvkit_hrd = CNVKIT_CALL.out.cnvkitsegment
            cnvkit_vcf_tumor = CNVKIT_EXPORT_VCF.out.cnvkit_vcf
                .filter { it -> it[2] == "full" && it[1].type == "T" }
                .map { val -> tuple(val[0], val[1], val[3]) }
        }
        else {
            CNVKIT_FULL ( bam_umi, params.cnvkit_reference, "full" )
            ch_versions = ch_versions.mix(CNVKIT_FULL.out.versions)

            CNVKIT_EXONS ( bam_umi, params.cnvkit_reference_exons, "exons" )
            ch_versions = ch_versions.mix(CNVKIT_EXONS.out.versions)

            CNVKIT_BACKBONE ( bam_umi, params.cnvkit_reference_backbone, "backbone" )
            ch_versions = ch_versions.mix(CNVKIT_BACKBONE.out.versions)

            cnvkit_backbone_input = CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3])
                .combine(germline_variants, by:[0])
                .map { group, meta, part, cns, cnr, germline_meta, vcf, tbi ->
                    tuple(group, meta, part, cns, cnr, vcf, tbi)
                }

            CNVKIT_PLOT ( cnvkit_backbone_input )
            ch_versions = ch_versions.mix(CNVKIT_PLOT.out.versions)

            cnvkit_split_input = CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3])
                .mix(
                    CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]),
                    CNVKIT_FULL.out.cnvkit_cns.join(CNVKIT_FULL.out.cnvkit_cnr, by:[0,1,3])
                )
                .combine(germline_variants, by:[0])
                .map { group, meta, part, cns, cnr, germline_meta, vcf, tbi ->
                    tuple(group, meta, part, cns, cnr, vcf, tbi)
                }

            cnvkit_split_call_input = cnvkit_split_input.map { group, meta, part, cns, cnr, vcf, tbi ->
                tuple(group, meta, part, cns, vcf, tbi)
            }

            cnvkit_nexus_input = cnvkit_split_input.map { group, meta, part, cns, cnr, vcf, tbi ->
                tuple(group, meta, part, cnr, vcf, tbi)
            }

            CNVKIT_CALL ( cnvkit_split_call_input, "false" )
            ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

            CNVKIT_EXPORT_VCF ( CNVKIT_CALL.out.cnvkitsegment )
            ch_versions = ch_versions.mix(CNVKIT_EXPORT_VCF.out.versions)

            CNVKIT_EXPORT_NEXUS_OGT ( cnvkit_nexus_input )
            ch_versions = ch_versions.mix(CNVKIT_EXPORT_NEXUS_OGT.out.versions)

            cnvkit_tc_input = CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3])
                .combine(germline_variants, by:[0])
                .map { group, meta, part, cns, cnr, germline_meta, vcf, tbi ->
                    tuple(group, meta, part, cns, cnr, vcf, tbi)
                }

            cnvkit_tc_call_input = cnvkit_tc_input.filter { group, meta, part, cns, cnr, vcf, tbi ->
                meta.purity && meta.purity.toString() != 'false'
            }.map { group, meta, part, cns, cnr, vcf, tbi ->
                tuple(group, meta, part, cns, vcf, tbi)
            }

            CNVKIT_CALL_TC ( cnvkit_tc_call_input, "true" )
            ch_versions = ch_versions.mix(CNVKIT_CALL_TC.out.versions)

            CNVKIT_EXPORT_VCF_TC ( CNVKIT_CALL_TC.out.cnvkitsegment )
            ch_versions = ch_versions.mix(CNVKIT_EXPORT_VCF_TC.out.versions)

            cnvkit_split_gens_input = CNVKIT_EXONS.out.cnvkit_cnr.mix(CNVKIT_BACKBONE.out.cnvkit_cnr)
                .combine(germline_variants, by:[0])
                .map { group, meta, cnr, part, germline_meta, vcf, tbi ->
                    tuple(group, meta, cnr, part, vcf, tbi)
                }

            CNVKIT_GENS ( cnvkit_split_gens_input )
            ch_versions = ch_versions.mix(CNVKIT_GENS.out.versions)

            ch_cnvkit_baf = CNVKIT_GENS.out.cnvkit_gens.map { group, meta, part, baf, cov ->
                tuple(group, meta, part, baf)
            }
            ch_cnvkit_cov = CNVKIT_GENS.out.cnvkit_gens.map { group, meta, part, baf, cov ->
                tuple(group, meta, part, cov)
            }
            CNVKIT_BAF_BGZIPTABIX ( ch_cnvkit_baf )
            CNVKIT_COV_BGZIPTABIX ( ch_cnvkit_cov )
            cnvkit_gens_compressed = CNVKIT_BAF_BGZIPTABIX.out.gz_index
                .join(CNVKIT_COV_BGZIPTABIX.out.gz_index, by:[0,1,2])
                .map { group, meta, part, baf_gz, baf_tbi, cov_gz, cov_tbi ->
                    tuple(group, meta, baf_gz, cov_gz)
                }

            MERGE_GENS ( cnvkit_gens_compressed.groupTuple(by:[0,1]) )
            ch_versions = ch_versions.mix(CNVKIT_BAF_BGZIPTABIX.out.versions)
            ch_versions = ch_versions.mix(CNVKIT_COV_BGZIPTABIX.out.versions)
            ch_versions = ch_versions.mix(MERGE_GENS.out.versions)

            cnvkitplot = CNVKIT_PLOT.out.cnvkitplot.filter { it -> it[2] == "backbone" }
            cnvkit_hrd = CNVKIT_CALL_TC.out.cnvkitsegment.mix(
                CNVKIT_CALL.out.cnvkitsegment.filter { group, meta, part, cns ->
                    part == "backbone" && (!meta.purity || meta.purity.toString() == 'false')
                }
            )
            cnvkit_vcf_tumor = CNVKIT_EXPORT_VCF.out.cnvkit_vcf.mix(CNVKIT_EXPORT_VCF_TC.out.cnvkit_vcf)
                .filter { it -> it[2] == "full" && it[1].type == "T" }
                .map { val -> tuple(val[0], val[1], val[3]) }
        }

    emit:
        plot      = cnvkitplot
        hrd       = cnvkit_hrd
        tumor_vcf = cnvkit_vcf_tumor
        gens      = MERGE_GENS.out.dbload
        gens_v4   = MERGE_GENS.out.gens_v4
        versions  = ch_versions
}
