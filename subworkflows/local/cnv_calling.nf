#!/usr/bin/env nextflow

include { CNVKIT                } from '../../modules/local/cnvkit/main'


workflow CNV_CALLING {
    take: 
        bam_umi
        concat_vcf

    main:
        CNVKIT ( bam_umi.combine(concat_vcf.filter { item -> item[1] == 'freebayes' }, by:[0]) )


    emit:
        baflogr = CNVKIT.out.baflogr

}