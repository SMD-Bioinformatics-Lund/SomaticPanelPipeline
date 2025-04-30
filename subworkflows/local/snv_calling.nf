#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { TNSCOPE                  } from '../../modules/local/sentieon/main'
include { PINDEL_CONFIG            } from '../../modules/local/pindel/main'
include { PINDEL_CALL              } from '../../modules/local/pindel/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { CONCATENATE_VCFS_BCFTOOLS} from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS as VT_AGG } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS as BT_AGG } from '../../modules/local/concatenate_vcfs/main'
include { MELT                     } from '../../modules/local/melt/main'
include { SVDB_MERGE_SINGLES       } from '../../modules/local/svdb/main'
include { BEDTOOLS_INTERSECT       } from '../../modules/local/filters/main'

workflow SNV_CALLING {
    take: 
        bam_umi                 // channel: [mandatory] [ val(group), val(meta), file("umi.bam"), file("umi.bam.bai"), file(bqsr) ]
        bam_dedup               // channel: [mandatory] [ val(group), val(meta), file(bam), file(bai)]
        beds                    // channel: [mandatory] [ file(bed) ]
        meta                    // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
        qc_values               // channel: [mandatory] [ val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) ]
        dedup_bam_is_metrics    // channel: [mandatory] [ val(group), val(meta), file(dedup_bam), file(dedup_bai), file(is_metrics) ]

    main:
        ch_versions = Channel.empty()

        // Pindel
        PINDEL_CONFIG ( dedup_bam_is_metrics )
        PINDEL_CALL ( dedup_bam_is_metrics, PINDEL_CONFIG.out.pindel_config )
        ch_versions         = ch_versions.mix(PINDEL_CALL.out.versions)

        // Variantcallers //
        // split by bed-file to speed up calling //
        FREEBAYES ( bam_umi, beds)
        ch_versions         = ch_versions.mix(FREEBAYES.out.versions.first())

        VARDICT ( bam_umi, beds)
        ch_versions         = ch_versions.mix(VARDICT.out.versions.first())

        TNSCOPE ( bam_umi, beds)
        ch_versions         = ch_versions.mix(TNSCOPE.out.versions.first())

        MELT ( bam_dedup.join(qc_values, by:[0,1])  )
        ch_versions         = ch_versions.mix(MELT.out.versions.first())
        BEDTOOLS_INTERSECT ( 
            MELT.out.melt_vcf,
            params.regions_bed
        )
        ch_versions         = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())

        if ( meta.filter( it -> it[1].type == "N" ) ) {
            SVDB_MERGE_SINGLES ( BEDTOOLS_INTERSECT.out.vcf_intersected.groupTuple() )
            MELT_MERGED     = SVDB_MERGE_SINGLES.out.singles_merged_vcf
            ch_versions     = ch_versions.mix(SVDB_MERGE_SINGLES.out.versions.first())
        }
        else {
            MELT_MERGED     = BEDTOOLS_INTERSECT.out.vcf_intersected
        }

        // Prepare vcf parts for concatenation //
        vcfparts_freebayes  = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        vcfparts_vardict    = VARDICT.out.vcfparts_vardict.groupTuple(by:[0,1])
        vcfparts_tnscope    = TNSCOPE.out.vcfparts_tnscope.groupTuple(by:[0,1])

        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat      = vcfparts_freebayes.mix(vcfparts_vardict,vcfparts_tnscope)

        // Join vcfs split by bedparts //
        CONCATENATE_VCFS { vcfs_to_concat }
        ch_versions         = ch_versions.mix(CONCATENATE_VCFS.out.versions.first())

        CONCATENATE_VCFS_BCFTOOLS { vcfs_to_concat }

        // Aggregate all callers to one VCF
        VT_AGG ( CONCATENATE_VCFS.out.concatenated_vcfs.mix(PINDEL_CALL.out.pindel_vcf,MELT_MERGED).groupTuple().join(meta.groupTuple()), "vt" )

        BT_AGG ( CONCATENATE_VCFS_BCFTOOLS.out.concatenated_vcfs.mix(PINDEL_CALL.out.pindel_vcf,MELT_MERGED).groupTuple().join(meta.groupTuple()), "bcftools" )

        ch_versions         = ch_versions.mix(VT_AGG.out.versions.first())

    emit:
        concat_vcfs =   CONCATENATE_VCFS.out.concatenated_vcfs  // channel: [ val(group), val(vc), file(vcf.gz) ]
        agg_vcf     =   BT_AGG.out.vcf_concat                   // channel: [ val(group), val(meta), file(agg.vcf) ]
        versions    =   ch_versions                             // channel: [ file(versions) ]

}