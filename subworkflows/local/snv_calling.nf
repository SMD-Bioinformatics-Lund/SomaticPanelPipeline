#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { TNSCOPE                  } from '../../modules/local/sentieon/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS           } from '../../modules/local/concatenate_vcfs/main'
include { MELT                     } from '../../modules/local/melt/main'
include { SVDB_MERGE_SINGLES       } from '../../modules/local/svdb/main'

workflow SNV_CALLING {
    take: 
        bam_umi             // tuple val(group), val(meta), file("umi.bam"), file("umi.bam.bai")
        bam_dedup           // tuple val(group), val(meta), file("dedup.bam"), file("dedup.bam.bai")
        beds
        meta
        qc_values           // tuple val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

    main:
        // Variantcallers //
        // split by bed-file to speed up calling //
        FREEBAYES ( bam_umi, beds)
        VARDICT ( bam_umi, beds)
        TNSCOPE ( bam_umi, beds)
        MELT ( bam_dedup.join(qc_values, by:[0,1])  )
        if ( meta.filter( it -> it[1].type == "N" ) ) {
            SVDB_MERGE_SINGLES ( MELT.out.melt_vcf.groupTuple() )
            MELT_MERGED = SVDB_MERGE_SINGLES.out.singles_merged_vcf
        }
        else {
            MELT_MERGED = MELT.out.melt_vcf
        }
        // Prepare vcf parts for concatenation //
        vcfparts_freebayes = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        vcfparts_vardict   = VARDICT.out.vcfparts_vardict.groupTuple(by:[0,1])
        vcfparts_tnscope   = TNSCOPE.out.vcfparts_tnscope.groupTuple(by:[0,1])
        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict,vcfparts_tnscope)
        // Join vcfs split by bedparts //
        CONCATENATE_VCFS { vcfs_to_concat }
        // Aggregate all callers to one VCF
        AGGREGATE_VCFS { CONCATENATE_VCFS.out.concatenated_vcfs.mix(MELT_MERGED).groupTuple().join(meta.groupTuple()) }

    emit:
        concat_vcfs = CONCATENATE_VCFS.out.concatenated_vcfs
        agg_vcf = AGGREGATE_VCFS.out.vcf_concat

}