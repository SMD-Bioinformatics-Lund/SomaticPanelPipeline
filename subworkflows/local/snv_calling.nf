#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { TNSCOPE                  } from '../../modules/local/sentieon/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS           } from '../../modules/local/concatenate_vcfs/main'
include { MELT                     } from '../../modules/local/melt/main'
include { SVDB_MERGE_SINGLES       } from '../../modules/local/svdb/main'
include { BEDTOOLS_INTERSECT       } from '../../modules/local/filters/main'

workflow SNV_CALLING {
    take: 
        bam_umi             // tuple val(group), val(meta), file("umi.bam"), file("umi.bam.bai")
        bam_dedup           // tuple val(group), val(meta), file("dedup.bam"), file("dedup.bam.bai")
        beds
        meta
        qc_values           // tuple val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

    main:
        ch_versions = Channel.empty()

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

        // Aggregate all callers to one VCF
        AGGREGATE_VCFS { CONCATENATE_VCFS.out.concatenated_vcfs.mix(MELT_MERGED).groupTuple().join(meta.groupTuple()) }
        ch_versions         = ch_versions.mix(AGGREGATE_VCFS.out.versions.first())

    emit:
        concat_vcfs =   CONCATENATE_VCFS.out.concatenated_vcfs
        agg_vcf     =   AGGREGATE_VCFS.out.vcf_concat
        versions    =   ch_versions

}