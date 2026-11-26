#!/usr/bin/env nextflow

include { PURECN_FORMAT_VCF                         } from '../../modules/local/purecn/main'
include { PURECN_COVERAGE                           } from '../../modules/local/purecn/main'
include { PURECN_CALL                               } from '../../modules/local/purecn/main'
include { PURECN_LOH_EXTRACT                        } from '../../modules/local/purecn/main'
include { PURECN_CNLOH_GENS                         } from '../../modules/local/purecn/main'

workflow PURECN {
    take:
        vcf_marked_germlines           // channel: [mandatory] [ val(group), val(meta), file(agg.pon.vep.markgerm.vcf) ]
        bam_umi                      // channel: [mandatory] [ val(group), val(meta), file(bam), file(bai) ]
        meta                         // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()
        ch_purecn_cnloh_snp = Channel.empty()
        ch_purecn_cnloh_gens = Channel.empty()
        ch_purecn_cnloh_baf_plot = Channel.empty()

        if (params.purecn) {
            PURECN_FORMAT_VCF ( vcf_marked_germlines )
            ch_versions = ch_versions.mix(PURECN_FORMAT_VCF.out.versions)

            tumor_bam_umi = bam_umi.filter { it[1].type == "T" }

            PURECN_COVERAGE ( tumor_bam_umi )
            ch_versions = ch_versions.mix(PURECN_COVERAGE.out.versions)

            PURECN_CALL ( 
                PURECN_FORMAT_VCF.out.vcf_purecn, 
                PURECN_COVERAGE.out.tumor_cov 
            )
            ch_versions = ch_versions.mix(PURECN_CALL.out.versions)

            PURECN_LOH_EXTRACT ( 
                PURECN_CALL.out.purecn_loh_raw_calls, 
                PURECN_CALL.out.purecn_log, 
                PURECN_CALL.out.purecn_rds
            )
            ch_versions = ch_versions.mix(PURECN_LOH_EXTRACT.out.versions)

            PURECN_CNLOH_GENS ( 
                PURECN_LOH_EXTRACT.out.purecn_loh_filtered_calls, 
                PURECN_CALL.out.purecn_snv_posteriors
            )
            ch_purecn_cnloh_snp = PURECN_CNLOH_GENS.out.purecn_cnloh_snp
            ch_purecn_cnloh_gens = PURECN_CNLOH_GENS.out.purecn_cnloh_gens
            ch_purecn_cnloh_baf_plot = PURECN_CNLOH_GENS.out.purecn_cnloh_baf_plot
            ch_versions = ch_versions.mix(PURECN_CNLOH_GENS.out.versions)
        }

    
    emit:
        purecn_cnloh_snp          = ch_purecn_cnloh_snp
        purecn_cnloh_gens         = ch_purecn_cnloh_gens
        purecn_cnloh_baf_plot     = ch_purecn_cnloh_baf_plot
        versions                  = ch_versions

}