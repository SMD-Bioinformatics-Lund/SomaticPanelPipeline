#!/usr/bin/env nextflow

include { PURECN_FORMAT_VCF                         } from '../../modules/local/purecn/main'
include { PURECN_COVERAGE                           } from '../../modules/local/purecn/main'
include { PURECN_CALL                               } from '../../modules/local/purecn/main'
include { PURECN_LOH_EXTRACT                        } from '../../modules/local/purecn/main'
include { PURECN_CNLOH_GENS                         } from '../../modules/local/purecn/main'

workflow PURECN {
    take:
        vcf_marked_germlines           // channel: [mandatory] [ val(group), val(meta), file(agg.pon.vep.markgerm.vcf) ]
        bam_umi                      // channel: [mandatory] [ val(group), val(meta), file(umi_bam), file(umi_bai), file(bqsr.table) ]

    main:
        ch_versions = Channel.empty()
        ch_purecn_cnloh_snp = Channel.empty()
        ch_purecn_cnloh_gens = Channel.empty()
        ch_purecn_cnloh_baf_plot = Channel.empty()

        if (params.purecn) {
            PURECN_FORMAT_VCF ( vcf_marked_germlines )
            ch_versions = ch_versions.mix(PURECN_FORMAT_VCF.out.versions)

            tumor_bam_umi = bam_umi
                .filter { it[1].type == "T" }
                .map { group, meta, bam, bai, bqsr ->
                    tuple(group, meta, bam, bai)}

            PURECN_COVERAGE ( tumor_bam_umi )
            ch_versions = ch_versions.mix(PURECN_COVERAGE.out.versions)

            ch_call_input = PURECN_FORMAT_VCF.out.vcf_purecn
                .map { group, meta, vcf, tbi -> tuple(group, vcf, tbi) }
                .join(PURECN_COVERAGE.out.tumor_cov)
                // -> [group, vcf, tbi, meta_tumor, cov]
            ch_call_input.view()
            
            PURECN_CALL ( ch_call_input)
            ch_versions = ch_versions.mix(PURECN_CALL.out.versions)

            ch_loh_extract_input = PURECN_CALL.out.purecn_loh_raw_calls
                .join(PURECN_CALL.out.purecn_log, by:[0,1])
                .join(PURECN_CALL.out.purecn_rds, by:[0,1])
                // -> [group, meta, loh_raw, log, rds]

            PURECN_LOH_EXTRACT ( ch_loh_extract_input )
            ch_versions = ch_versions.mix(PURECN_LOH_EXTRACT.out.versions)

            ch_gens_input = PURECN_LOH_EXTRACT.out.purecn_loh_filtered
                .join(PURECN_CALL.out.purecn_snv_posteriors, by:[0,1])
                // -> [group, meta, loh_filtered, posteriors]

            PURECN_CNLOH_GENS ( ch_gens_input)
            ch_versions = ch_versions.mix(PURECN_CNLOH_GENS.out.versions)

            ch_purecn_cnloh_snp = PURECN_CNLOH_GENS.out.purecn_cnloh_snp
            ch_purecn_cnloh_gens = PURECN_CNLOH_GENS.out.purecn_cnloh_gens
            ch_purecn_cnloh_baf_plot = PURECN_CNLOH_GENS.out.purecn_cnloh_baf_plot

        }

    
    emit:
        purecn_cnloh_snp          = ch_purecn_cnloh_snp
        purecn_cnloh_gens         = ch_purecn_cnloh_gens
        purecn_cnloh_baf_plot     = ch_purecn_cnloh_baf_plot
        versions                  = ch_versions

}