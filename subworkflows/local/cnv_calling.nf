#!/usr/bin/env nextflow

include { CNVKIT                } from '../../modules/local/cnvkit/main'
include { GATKCOV_BAF           } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT         } from '../../modules/local/GATK/main'
include { GATKCOV_CALL          } from '../../modules/local/GATK/main'


workflow CNV_CALLING {
    take: 
        bam_umi
        germline_variants
        meta
        concat_vcf

    main:
        CNVKIT { bam_umi.join(meta, by:[0,1,2]).join(germline_variants) }
        GATKCOV_BAF { bam_umi }
        GATKCOV_COUNT { bam_umi }
        GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count, by:[0,1,2]) }
        

//CNVKIT { bam_umi.join(meta, by:[0,1,2]).join(concat_vcf.filter { item -> item[1] == 'freebayes' } ).view() }
    emit:
        baflogr = CNVKIT.out.baflogr
        cnvkitsegment = CNVKIT.out.cnvkitsegment
        cnvkitsegment_purity = CNVKIT.out.cnvkitsegment_purity

}