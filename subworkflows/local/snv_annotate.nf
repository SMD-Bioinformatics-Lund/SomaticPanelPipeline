#!/usr/bin/env nextflow

include { PON_FILTER               } from '../../modules/local/filters/main'
include { FFPE_PON_FILTER          } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { MARK_GERMLINES           } from '../../modules/local/filters/main'
include { FILTER_FOR_CNV           } from '../../modules/local/filters/main'
include { CONTAMINATION            } from '../../modules/local/filters/main'
include { VCFANNO                  } from '../../modules/local/filters/main'

workflow SNV_ANNOTATE {
    take: 
        agg_vcf
        concat_vcfs
        meta

    main:
        // Filter with PoN, annotate with VEP, mark germlines
        PON_FILTER { agg_vcf }

        // OPTIONAL PARTS SET BY CONFIG T/F
            // FFPE-PON if set to true
            if (params.ffpe_pon) {
                FFPE_PON_FILTER { PON_FILTER.out.vcf_pon } // needs to be merged with normal PON above, myeloid should not be be FFPE annotated
                VCFANNO_INPUT = FFPE_PON_FILTER.out.vcf_pon_ffpe
            }
            else {
                VCFANNO_INPUT = PON_FILTER.out.vcf_pon
            }
            // ENIGMA if set to true
            if (params.vcfanno) {
                VCFANNO { VCFANNO_INPUT }
                PON_VEP = VCFANNO.out.vcf_enigma
            }
            else {
                PON_VEP = VCFANNO_INPUT
            }
        /////////////////////////////////////////

        // NON-OPTIONAL, NEEDED BY COYOTE
        ANNOTATE_VEP { PON_VEP } 
        MARK_GERMLINES { ANNOTATE_VEP.out.vcf_vep }
        // BAF for CNVkit //
        FILTER_FOR_CNV { ANNOTATE_VEP.out.vcf_vep.join(concat_vcfs.filter { item -> item[1] == 'freebayes' })  }
        // contamination values from VCF //
        CONTAMINATION { ANNOTATE_VEP.out.vcf_vep }

    emit:
        germline_variants = FILTER_FOR_CNV.out.vcf_only_germline
        finished_vcf = MARK_GERMLINES.out.vcf_germline

}