#!/usr/bin/env nextflow

include { CNVKIT                } from '../../modules/local/cnvkit/main'


workflow CNV_CALLING {
    take: 
        bam_umi
        concat_vcf
        meta

    main:
        CNVKIT ( bam_umi.join(meta, by:[0,1,2]).combine(concat_vcf.filter { item -> item[1] == 'freebayes' }, by:[0]) )


    emit:
        baflogr = CNVKIT.out.baflogr
        cnvkitsegment = CNVKIT.out.cnvkitsegment
        cnvkitsegment_purity = CNVKIT.out.cnvkitsegment_purity

}