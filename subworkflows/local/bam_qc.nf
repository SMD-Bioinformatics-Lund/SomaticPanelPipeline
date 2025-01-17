#!/usr/bin/env nextflow

include { SENTIEON_QC          } from '../../modules/local/sentieon/main'
include { QC_TO_CDM            } from '../../modules/local/qc/main'
include { LOWCOV               } from '../../modules/local/qc/main'
include { QC_VALUES            } from '../../modules/local/qc/main'
include { VERIFYBAMID          } from '../../modules/local/verifybamid/main'
include { ALLELE_CALL          } from '../../modules/local/idSnp/main'
include { SNP_CHECK            } from '../../modules/local/idSnp/main'
include { PAIRGEN_CDM          } from '../../modules/local/idSnp/main'
include { LOWCOV_D4            } from '../../modules/local/qc/main'

workflow BAM_QC {
    take:        
        bam_umi         // channel: [mandatory] [ val(group), val(meta), file("umi.bam"), file("umi.bam.bai"), file(bqsr) ]
        bam_dedup       // channel: [ val(group), val(meta), file(bam), file(bai) ]
        dedup_metrics   // channel: [ val(group), val(meta), file(dedup_metrics) ]

    main:
        ch_versions = Channel.empty()
        
        SENTIEON_QC ( bam_dedup.join(dedup_metrics, by:[0,1]) )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)
        
        QC_VALUES ( SENTIEON_QC.out.qc_cdm )

        QC_TO_CDM ( SENTIEON_QC.out.qc_cdm )

        LOWCOV ( bam_dedup )
        ch_versions = ch_versions.mix(LOWCOV.out.versions)

        LOWCOV_D4 ( bam_dedup )
        ch_versions = ch_versions.mix(LOWCOV_D4.out.versions)

        // Check genotypes of ID-SNPs
        ALLELE_CALL (bam_dedup)
        ch_versions = ch_versions.mix(ALLELE_CALL.out.versions)

        SNP_CHECK(ALLELE_CALL.out.sample_id_genotypes.groupTuple())
        ch_versions = ch_versions.mix(SNP_CHECK.out.versions)

        PAIRGEN_CDM (SNP_CHECK.out.idsnp_checked)

        // // Calculate cross-sample contamination
        // VERIFYBAMID { bam_umi }
        // ch_versions = ch_versions.mix(VERIFYBAMID.out.versions)
    emit:
        qcdone                  =   QC_TO_CDM.out.cdm_done                  // channel: [ val(group), val(meta), file(cdm) ]
        lowcov                  =   LOWCOV.out.lowcov_regions               // channel: [ val(group), val(meta.type), file(lowcov.bed) ]
        lowcov_d4               =   LOWCOV_D4.out.coyote_cov_json           // channel: [ val(group), file(cov.json) ]
        melt_qc                 =   QC_VALUES.out.qc_melt_val               // channel: [ val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) ]
        versions                =   ch_versions                             // channel: [ file(versions) ]
        dedup_bam_is_metrics    =   SENTIEON_QC.out.dedup_bam_is_metrics    // channel: [ val(group), val(meta), file(is_metrics.txt) ] 
}