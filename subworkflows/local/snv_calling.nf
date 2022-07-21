#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS           } from '../../modules/local/concatenate_vcfs/main'

workflow SNV_CALLING {
    take: 
        bam_umi
        beds
        meta

    main:
        // Variantcallers //
        FREEBAYES ( bam_umi, beds)
        VARDICT ( bam_umi, beds)
        // Prepare vcf parts for concatenation
        vcfparts_freebayes = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        vcfparts_vardict   = VARDICT.out.vcfparts_vardict.groupTuple(by:[0,1])
        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict)
        // Join vcfs split by bedparts
        CONCATENATE_VCFS { vcfs_to_concat }
        // Aggregate all callers to one VCF
        AGGREGATE_VCFS { CONCATENATE_VCFS.out.concatenated_vcfs.join(meta).groupTuple() }

    emit:
        concat_vcfs = CONCATENATE_VCFS.out.concatenated_vcfs

}