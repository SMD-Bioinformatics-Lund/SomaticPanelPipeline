#!/usr/bin/env nextflow

include { FREEBAYES                             } from '../../modules/local/freebayes/main'
include { FILTER as FREEBAYES_FILTER_LOWCOV    } from '../../modules/local/vcflib/main'
include { FILTER as FREEBAYES_FILTER_LOWFRQ    } from '../../modules/local/vcflib/main'
include { GLXGT as FREEBAYES_VCFGLXGT          } from '../../modules/local/vcflib/main'
include { FILTER_FREEBAYES                     } from '../../modules/local/filters/main'
include { ANNOTATE as FREEBAYES_REMOVE_AD      } from '../../modules/local/bcftools/main'
include { VARDICT                               } from '../../modules/local/vardict/main'
include { FILTER_VARDICT                        } from '../../modules/local/filters/main'
include { TNSCOPE                               } from '../../modules/local/sentieon/main'
include { PINDEL_CONFIG                         } from '../../modules/local/pindel/main'
include { PINDEL_CALL                           } from '../../modules/local/pindel/main'
include { PINDEL_TO_VCF                         } from '../../modules/local/pindel/main'
include { FILTER_PINDEL                         } from '../../modules/local/filters/main'
include { VCF_CONCAT                            } from '../../modules/local/vcftools/main'
include { VCF_SORT as VCF_CONCATENATED_SORT     } from '../../modules/local/vcftools/main'
include { TABIX_BGZIPTABIX                      } from '../../modules/local/tabix/main'
include { BCFTOOLS_NORM                         } from '../../modules/local/bcftools/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_EXACT  } from '../../modules/local/bcftools/main'
include { CONCATENATE_VCFS                      } from '../../modules/local/concatenate_vcfs/main'
// include { CONCATENATE_VCFS_BCFTOOLS             } from '../../modules/local/concatenate_vcfs/main'
// include { AGGREGATE_VCFS as VT_AGG              } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS as BT_AGG              } from '../../modules/local/concatenate_vcfs/main'
include { VCF_SORT as VCF_AGGREGATE_SORT        } from '../../modules/local/vcftools/main'
include { MELT                                  } from '../../modules/local/melt/main'
include { SVDB_MERGE_SINGLES                    } from '../../modules/local/svdb/main'
include { BEDTOOLS_INTERSECT                    } from '../../modules/local/bedtools/main'
include { FILTER_TNSCOPE                        } from '../../modules/local/filters/main'

workflow SNV_CALLING {
    take:
        bam_umi                 // channel: [mandatory] [ val(group), val(meta), file("umi.bam"), file("umi.bam.bai"), file(bqsr) ]
        bam_dedup               // channel: [mandatory] [ val(group), val(meta), file(bam), file(bai)]
        beds                    // channel: [mandatory] [ file(bed) ]
        // meta                    // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        qc_values               // channel: [mandatory] [ val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) ]
        dedup_bam_is_metrics    // channel: [mandatory] [ val(group), val(meta), file(dedup_bam), file(dedup_bai), file(is_metrics) ]

    main:
        ch_versions = Channel.empty()

        // Pindel
        PINDEL_CONFIG ( dedup_bam_is_metrics )
        PINDEL_CALL ( dedup_bam_is_metrics, PINDEL_CONFIG.out.pindel_config )
        ch_versions         = ch_versions.mix(PINDEL_CALL.out.versions)

        PINDEL_TO_VCF ( PINDEL_CALL.out.pindel_raw )
        ch_versions         = ch_versions.mix(PINDEL_TO_VCF.out.versions)

        FILTER_PINDEL ( PINDEL_TO_VCF.out.pindel_vcf_unfiltered )
        ch_versions         = ch_versions.mix(FILTER_PINDEL.out.versions.first())

        // Variantcallers //
        // split by bed-file to speed up calling //
        FREEBAYES ( bam_umi, beds)
        FREEBAYES_FILTER_LOWCOV ( FREEBAYES.out.vcfparts_freebayes )
        FREEBAYES_FILTER_LOWFRQ ( FREEBAYES_FILTER_LOWCOV.out.filtered_vcf )
        FREEBAYES_VCFGLXGT ( FREEBAYES_FILTER_LOWFRQ.out.filtered_vcf )
        FILTER_FREEBAYES ( FREEBAYES_VCFGLXGT.out.fixed_genotypes_vcf )
        FREEBAYES_REMOVE_AD ( FILTER_FREEBAYES.out.filtered_vcf )
        ch_versions         = ch_versions.mix(FREEBAYES.out.versions.first())
        ch_versions         = ch_versions.mix(FREEBAYES_FILTER_LOWCOV.out.versions.first())
        ch_versions         = ch_versions.mix(FREEBAYES_FILTER_LOWFRQ.out.versions.first())
        ch_versions         = ch_versions.mix(FREEBAYES_VCFGLXGT.out.versions.first())
        ch_versions         = ch_versions.mix(FILTER_FREEBAYES.out.versions.first())
        ch_versions         = ch_versions.mix(FREEBAYES_REMOVE_AD.out.versions.first())

        VARDICT ( bam_umi, beds)
        FILTER_VARDICT ( VARDICT.out.raw_vcfparts_vardict )
        ch_versions         = ch_versions.mix(VARDICT.out.versions.first())
        ch_versions         = ch_versions.mix(FILTER_VARDICT.out.versions.first())

        TNSCOPE ( bam_umi, beds)
        FILTER_TNSCOPE ( TNSCOPE.out.vcfparts_tnscope )
        ch_versions         = ch_versions.mix(TNSCOPE.out.versions.first())


        //TODO: THis is not run by any profile, should we remove this?
        MELT ( bam_dedup.join(qc_values, by:[0,1])  )
        ch_versions         = ch_versions.mix(MELT.out.versions.first())

        //TODO: THis is not run by any profile, should we remove this?
        ch_melt_intersect_input = MELT.out.melt_vcf.map { group, meta, vc, vcf ->
            tuple(group, meta, vc, vcf)
        }
        ch_regions_bed = Channel.value(file(params.regions_bed))

        BEDTOOLS_INTERSECT (
            ch_melt_intersect_input,
            ch_regions_bed
        )
        ch_versions         = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())


        /////////////////////////////////////////////////////////////////////////////////
        ch_melt_grouped = BEDTOOLS_INTERSECT.out.vcf_intersected
            .groupTuple(by:[0,2])
            .map { group, metas, vc, vcfs ->
                [
                    group,
                    metas,
                    vc,
                    vcfs.sort { a, b -> a.name <=> b.name }
                ]
            }

        ch_melt_split = ch_melt_grouped.branch { group, metas, vc, vcfs ->
            has_normal:
                metas.any { sample_meta -> sample_meta.type == 'N' || sample_meta.type == 'normal' }

            tumor_only:
                true
        }

        SVDB_MERGE_SINGLES ( ch_melt_split.has_normal )

        MELT_MERGED = SVDB_MERGE_SINGLES.out.singles_merged_vcf.mix(
            ch_melt_split.tumor_only.map { group, metas, vc, vcfs ->
                [
                    group,
                    metas,
                    vc,
                    vcfs[0]
                ]
            }
        )

        //////////////////////////////////////////////////////////////////////////////////

        // if ( meta.filter( it -> it[1].type == "N" ) ) {
        //     SVDB_MERGE_SINGLES ( BEDTOOLS_INTERSECT.out.vcf_intersected.groupTuple() )
        //     MELT_MERGED     = SVDB_MERGE_SINGLES.out.singles_merged_vcf
        //     ch_versions     = ch_versions.mix(SVDB_MERGE_SINGLES.out.versions.first())
        // }
        // else {
        //     MELT_MERGED     = BEDTOOLS_INTERSECT.out.vcf_intersected
        // }

        // Prepare vcf parts for concatenation //

        vcfparts_freebayes  = FREEBAYES_REMOVE_AD.out.annotated_vcf.groupTuple(by:[0,2]).map { group, metas, vc, vcfparts ->
            [
                group,
                metas[0], // all metas should be the same for a given group and vc, so just take the first one
                vc,
                vcfparts.sort{ a, b -> a.name <=> b.name } // sort by filename to ensure same order for all callers
            ]
        }
        vcfparts_vardict    = FILTER_VARDICT.out.vcfparts_vardict.groupTuple(by:[0,2]).map { group, metas, vc, vcfparts ->
            [
                group,
                metas[0], // all metas should be the same for a given group and vc, so just take the first one
                vc,
                vcfparts.sort{ a, b -> a.name <=> b.name } // sort by filename to ensure same order for all callers
            ]
        }
        vcfparts_tnscope    = FILTER_TNSCOPE.out.vcfparts_tnscope_filtered.groupTuple(by:[0,2]).map { group, metas, vc, vcfparts ->
            [
                group,
                metas[0], // all metas should be the same for a given group and vc, so just take the first one
                vc,
                vcfparts.sort{ a, b -> a.name <=> b.name } // sort by filename to ensure same order for all callers
            ]
        }

        // vcfparts_freebayes.view()
        // vcfparts_vardict.view()
        // vcfparts_tnscope.view()

        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat      = vcfparts_freebayes.mix(vcfparts_vardict,vcfparts_tnscope)

        // concatenate vcfs by bedparts
        VCF_CONCAT ( vcfs_to_concat )
        ch_versions         = ch_versions.mix(VCF_CONCAT.out.versions.first())


        VCF_CONCATENATED_SORT ( VCF_CONCAT.out.concatenated_vcf )
        ch_versions         = ch_versions.mix(VCF_CONCATENATED_SORT.out.versions.first())


        TABIX_BGZIPTABIX( VCF_CONCATENATED_SORT.out.sorted_vcf )
        ch_versions         = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

        BCFTOOLS_NORM ( TABIX_BGZIPTABIX.out.gz_index )
        BCFTOOLS_NORM_EXACT ( BCFTOOLS_NORM.out.normalized_vcfs )


        // Join vcfs split by bedparts //
        // CONCATENATE_VCFS { vcfs_to_concat }
        // ch_versions         = ch_versions.mix(CONCATENATE_VCFS.out.versions.first())

        // CONCATENATE_VCFS_BCFTOOLS { vcfs_to_concat }

        // Aggregate all callers to one VCF
        // VT_AGG ( CONCATENATE_VCFS.out.concatenated_vcfs.mix(PINDEL_CALL.out.pindel_vcf,MELT_MERGED).groupTuple().join(meta.groupTuple()), "vt" )

        BT_AGG (
            BCFTOOLS_NORM_EXACT.out.normalized_vcfs.map {
                group, meta, vc, vcf, tbi ->
                [
                    group,
                    meta,
                    vc,
                    vcf
                ]
            } .mix(
                FILTER_PINDEL.out.pindel_vcf,
                MELT_MERGED
            ).groupTuple().map {
                group, meta, vc, vcfs ->
                [
                    group,
                    meta[0],
                    vcfs.flatten()
                ]
            }
        )

        // BT_AGG ( BCFTOOLS_NORM_EXACT.out.normalized_vcfs.mix(PINDEL_GZ.out.gz_index,MELT_MERGED).groupTuple().join(meta.groupTuple()) )

        VCF_AGGREGATE_SORT (
            BT_AGG.out.vcf_agg.map {
                group, meta, vcf ->
                [
                    group,
                    meta,
                    [],
                    vcf
                ]
            }

        )

        ch_agg_vcf_sorted_meta = VCF_AGGREGATE_SORT.out.sorted_vcf.map {
            group, meta, vc, vcf ->
            [
                group,
                meta,
                vcf
            ]
        }

        ch_versions         = ch_versions.mix(BT_AGG.out.versions.first())

    emit:
        concat_vcfs =   BCFTOOLS_NORM_EXACT.out.normalized_vcfs             // channel: [ val(group), val(meta), val(vc), file(vcf.gz), file(vcf.gz.tbi) ]
        agg_vcf     =   ch_agg_vcf_sorted_meta                              // channel: [ val(group), val(meta), file(agg.vcf) ]
        versions    =   ch_versions                                         // channel: [ file(versions) ]

}
