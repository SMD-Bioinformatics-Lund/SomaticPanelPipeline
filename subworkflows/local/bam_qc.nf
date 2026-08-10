#!/usr/bin/env nextflow

include { SENTIEON_QC                               } from '../../modules/local/sentieon/main'
include { DEPTH as SAMBAMBA_PANEL_DEPTH             } from '../../modules/local/sambamba/main'
include { BEDTOOLS_INTERSECT as LOWCOV_INTERSECT    } from '../../modules/local/bedtools/main'
include { QC_TO_CDM                                 } from '../../modules/local/qc/main'
include { SENTIEON_QC_TO_CDM                        } from '../../modules/local/qc/main'
include { PANEL_DEPTH                               } from '../../modules/local/qc/main'
include { LOWCOV                                    } from '../../modules/local/qc/main'
include { QC_VALUES                                 } from '../../modules/local/qc/main'
include { VERIFYBAMID                               } from '../../modules/local/verifybamid/main'
include { ANNOTATE as BCFTOOLS_ANNOTATE_IDSNP       } from '../../modules/local/bcftools/main'
include { BCFTOOLS_MPILEUP_CALL                     } from '../../modules/local/bcftools/main'
include { BCFTOOLS_QUERY_IDSNP                      } from '../../modules/local/bcftools/main'
include { GENOTYPE2JSON                             } from '../../modules/local/idSnp/main'
include { SNP_CHECK                                 } from '../../modules/local/idSnp/main'
include { PAIRGEN_CDM                               } from '../../modules/local/idSnp/main'
include { LOWCOV_D4                                 } from '../../modules/local/qc/main'

workflow BAM_QC {
    take:        
        bam_umi         // channel: [mandatory] [ val(group), val(meta), file("umi.bam"), file("umi.bam.bai"), file(bqsr) ]
        bam_dedup       // channel: [ val(group), val(meta), file(bam), file(bai) ]
        dedup_metrics   // channel: [ val(group), val(meta), file(dedup_metrics) ]

    main:
        ch_versions = Channel.empty()
        
        SENTIEON_QC ( bam_dedup.join(dedup_metrics, by:[0,1]) )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)
        SENTIEON_QC_TO_CDM( SENTIEON_QC.out.qc_files.join(dedup_metrics, by:[0,1]) )
        
        QC_VALUES ( SENTIEON_QC_TO_CDM.out.qc_cdm )

        QC_TO_CDM ( SENTIEON_QC_TO_CDM.out.qc_cdm )

        LOWCOV_D4 ( bam_dedup )
        ch_versions = ch_versions.mix(LOWCOV_D4.out.versions)

        SAMBAMBA_PANEL_DEPTH ( bam_dedup )
        ch_versions = ch_versions.mix(SAMBAMBA_PANEL_DEPTH.out.versions)
        
        PANEL_DEPTH ( SAMBAMBA_PANEL_DEPTH.out.depth_tsv )
        ch_versions = ch_versions.mix(PANEL_DEPTH.out.versions)


        ch_lowcov_intersect_input = PANEL_DEPTH.out.lowcov_bed.map { group, meta, bed ->
            tuple(group, meta, 'lowcov', bed)
        }
        ch_gene_regions = Channel.value(file(params.gene_regions))

        LOWCOV_INTERSECT (
            ch_lowcov_intersect_input,
            ch_gene_regions
        )
        ch_versions = ch_versions.mix(LOWCOV_INTERSECT.out.versions)

        LOWCOV(LOWCOV_INTERSECT.out.intersected)
        ch_versions = ch_versions.mix(LOWCOV.out.versions)

        // Check genotypes of ID-SNPs
        BCFTOOLS_MPILEUP_CALL (bam_dedup)
        ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP_CALL.out.versions)

        BCFTOOLS_ANNOTATE_IDSNP (
            BCFTOOLS_MPILEUP_CALL.out.raw_vcf.map { group, meta, vcf ->
                tuple(group, meta, 'idsnp', vcf)
            }
        )
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_IDSNP.out.versions)

        BCFTOOLS_QUERY_IDSNP (
            BCFTOOLS_ANNOTATE_IDSNP.out.annotated_vcf.map { group, meta, vc, vcf ->
                tuple(group, meta, vcf)
            }
        )
        ch_versions = ch_versions.mix(BCFTOOLS_QUERY_IDSNP.out.versions)

        GENOTYPE2JSON (BCFTOOLS_QUERY_IDSNP.out.genotypes)
        ch_versions = ch_versions.mix(GENOTYPE2JSON.out.versions)

        SNP_CHECK(GENOTYPE2JSON.out.sample_id_genotypes.groupTuple())
        ch_versions = ch_versions.mix(SNP_CHECK.out.versions)

        PAIRGEN_CDM (SNP_CHECK.out.idsnp_tumor.mix(SNP_CHECK.out.idsnp_normal))

        // // Calculate cross-sample contamination
        // VERIFYBAMID { bam_umi }
        // ch_versions = ch_versions.mix(VERIFYBAMID.out.versions)
    emit:
        qcdone                  =   QC_TO_CDM.out.cdm_done                  // channel: [ val(group), val(meta), file(cdm) ]
        lowcov                  =   LOWCOV.out.lowcov_regions               // channel: [ val(group), val(meta.type), file(lowcov.bed/cov.json) ]
        lowcov_d4               =   LOWCOV_D4.out.coyote_cov_json
        melt_qc                 =   QC_VALUES.out.qc_melt_val               // channel: [ val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) ]
        versions                =   ch_versions                             // channel: [ file(versions) ]
        dedup_bam_is_metrics    =   SENTIEON_QC.out.dedup_bam_is_metrics    // channel: [ val(group), val(meta), file(is_metrics.txt) ] 
}
