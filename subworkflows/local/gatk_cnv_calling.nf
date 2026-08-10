#!/usr/bin/env nextflow

include { GATKCOV_BAF                                               } from '../../modules/local/GATK/main'
include { GATK_COLLECT_READ_COUNTS                                  } from '../../modules/local/GATK/main'
include { GATK_DENOISE_READ_COUNTS                                  } from '../../modules/local/GATK/main'
include { GATK_PLOT_DENOISED_COPY_RATIOS                            } from '../../modules/local/GATK/main'
include { GATK_MODEL_SEGMENTS                                       } from '../../modules/local/GATK/main'
include { GATK_CALL_COPY_RATIO_SEGMENTS                             } from '../../modules/local/GATK/main'
include { GATK_PLOT_MODELED_SEGMENTS                                } from '../../modules/local/GATK/main'
include { GATK2VCF                                                  } from '../../modules/local/GATK/main'
include { GATK_COUNT_GERMLINE                                       } from '../../modules/local/GATK/main'
include { GATK_CALL_PLOIDY                                          } from '../../modules/local/GATK/main'
include { GATK_CALL_GERMLINE_CNV                                    } from '../../modules/local/GATK/main'
include { GATK_EXTRACT_GERMLINE_CNV_SHARDS                          } from '../../modules/local/GATK/main'
include { FILTER_GATK                                               } from '../../modules/local/GATK/main'
include { MERGE_GATK                                                } from '../../modules/local/GATK/main'
include { MERGE_GATK_TUMOR                                          } from '../../modules/local/GATK/main'
include { POSTPROCESS                                               } from '../../modules/local/GATK/main'


workflow GATK {
    take:
        bam_umi      // channel: [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]
        gatk_ref     // channel: [ val(integer), val(part_of_genome) ]

    main:
        ch_versions = Channel.empty()

        GATKCOV_BAF ( bam_umi )
        ch_versions = ch_versions.mix(GATKCOV_BAF.out.versions)

        GATK_COLLECT_READ_COUNTS ( bam_umi )
        ch_versions = ch_versions.mix(GATK_COLLECT_READ_COUNTS.out.versions)

        GATK_DENOISE_READ_COUNTS ( GATK_COLLECT_READ_COUNTS.out.read_counts )
        ch_versions = ch_versions.mix(GATK_DENOISE_READ_COUNTS.out.versions)

        GATK_PLOT_DENOISED_COPY_RATIOS ( GATK_DENOISE_READ_COUNTS.out.denoised_counts )
        ch_versions = ch_versions.mix(GATK_PLOT_DENOISED_COPY_RATIOS.out.versions)

        GATK_MODEL_SEGMENTS { GATKCOV_BAF.out.gatk_baf.join(GATK_DENOISE_READ_COUNTS.out.denoised_counts,by:[0,1]).groupTuple() }
        ch_versions = ch_versions.mix(GATK_MODEL_SEGMENTS.out.versions)

        ch_model_segments = GATK_MODEL_SEGMENTS.out.model_segments.map { group, meta, cr_seg, hets, model_final ->
            tumor_meta = meta.find { it.type == "T" || it.type == "tumor" }
            tuple(group, tumor_meta, cr_seg, hets, model_final)
        }

        GATK_CALL_COPY_RATIO_SEGMENTS ( ch_model_segments )
        ch_versions = ch_versions.mix(GATK_CALL_COPY_RATIO_SEGMENTS.out.versions)

        ch_tumor_denoised = GATK_DENOISE_READ_COUNTS.out.denoised_counts
            .filter { group, meta, stdCR, denoised -> meta.type == "T" || meta.type == "tumor" }
            .map { group, meta, stdCR, denoised -> tuple(group, meta, denoised) }

        ch_modeled_plot_input = ch_tumor_denoised
            .join(ch_model_segments, by:[0,1])
            .map { group, meta, denoised, cr_seg, hets, model_final ->
                tuple(group, meta, denoised, cr_seg, hets, model_final)
            }

        GATK_PLOT_MODELED_SEGMENTS ( ch_modeled_plot_input )
        ch_versions = ch_versions.mix(GATK_PLOT_MODELED_SEGMENTS.out.versions)

        GATK2VCF ( GATK_CALL_COPY_RATIO_SEGMENTS.out.called_segments )
        ch_versions = ch_versions.mix(GATK2VCF.out.versions)

        MERGE_GATK_TUMOR ( GATK2VCF.out.tumor_vcf )
        ch_versions = ch_versions.mix(MERGE_GATK_TUMOR.out.versions)

        GATK_COUNT_GERMLINE ( bam_umi.filter { it -> it[1].type == "N" || it[1].type == "normal" })
        ch_versions = ch_versions.mix(GATK_COUNT_GERMLINE.out.versions)

        GATK_CALL_PLOIDY ( GATK_COUNT_GERMLINE.out.count_germline )
        ch_versions = ch_versions.mix(GATK_CALL_PLOIDY.out.versions)

        GATK_CALL_GERMLINE_CNV ( GATK_COUNT_GERMLINE.out.count_germline.join(GATK_CALL_PLOIDY.out.gatk_ploidy,by:[0,1]).groupTuple(by:[0,1]).combine(gatk_ref) )
        ch_versions = ch_versions.mix(GATK_CALL_GERMLINE_CNV.out.versions)

        CALLED = GATK_CALL_GERMLINE_CNV.out.gatk_call_germline.groupTuple(by:[0,1])
        PLOIDY = GATK_CALL_PLOIDY.out.gatk_ploidy.groupTuple(by:[0,1])

        GATK_EXTRACT_GERMLINE_CNV_SHARDS ( CALLED.join(PLOIDY,by:[0,1]).combine(gatk_ref.groupTuple(by:[3])) )
        ch_versions = ch_versions.mix(GATK_EXTRACT_GERMLINE_CNV_SHARDS.out.versions)

        POSTPROCESS ( GATK_EXTRACT_GERMLINE_CNV_SHARDS.out.extracted_germline_cnv )
        ch_versions = ch_versions.mix(POSTPROCESS.out.versions)

        FILTER_GATK ( POSTPROCESS.out.gatk_germline_segmentsvcf )
        ch_versions = ch_versions.mix(FILTER_GATK.out.versions)

        MERGE_GATK ( FILTER_GATK.out.gatk_filtered_vcf )
        ch_versions = ch_versions.mix(MERGE_GATK.out.versions)

    emit:
        gatcov_plot = GATK_PLOT_MODELED_SEGMENTS.out.modeled_plot
        tumor_vcf   = MERGE_GATK_TUMOR.out.tumor_vcf_merged
        normal_vcf  = MERGE_GATK.out.gatk_normal_vcf
        versions    = ch_versions
}
