#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'

workflow SNV_CALLING {
    take: 
        bam_umi
        beds

    main:
        FREEBAYES ( bam_umi, beds)
        // Prepare vcf parts for concatenation
        vcfparts_freebayes = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        //vcfparts_tnscope   = vcfparts_tnscope.groupTuple(by:[0,1])
        //vcfparts_vardict   = vcfparts_vardict.groupTuple(by:[0,1])
        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat = vcfparts_freebayes

        CONCATENATE_VCFS { vcfparts_freebayes }

    emit:
        concat_vcfs = CONCATENATE_VCFS.out.concatenated_vcfs

}