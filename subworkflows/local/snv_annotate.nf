#!/usr/bin/env nextflow

include { PON_FILTER               } from '../../modules/local/filters/main'
include { FFPE_PON_FILTER          } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { MARK_GERMLINES           } from '../../modules/local/filters/main'
include { FILTER_FOR_CNV           } from '../../modules/local/filters/main'
include { VCFANNO                  } from '../../modules/local/filters/main'

workflow SNV_ANNOTATE {
    take: 
        agg_vcf         // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        concat_vcfs     // channel: [mandatory] [ val(group), val(vc), file(vcf.gz) ]
        meta            // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        // Filter with PoN, annotate with VEP, mark germlines
        PON_FILTER { agg_vcf }
        ch_versions = ch_versions.mix(PON_FILTER.out.versions)

        // OPTIONAL PARTS SET BY CONFIG T/F
            // FFPE-PON if set to true
            if (params.ffpe_pon) {
                FFPE_PON_FILTER { PON_FILTER.out.vcf_pon } // needs to be merged with normal PON above, myeloid should not be be FFPE annotated
                VCFANNO_INPUT = FFPE_PON_FILTER.out.vcf_pon_ffpe
                ch_versions = ch_versions.mix(FFPE_PON_FILTER.out.versions)
            }
            else {
                VCFANNO_INPUT = PON_FILTER.out.vcf_pon
            }

            // ENIGMA if set to true
            if (params.vcfanno) {
                VCFANNO { VCFANNO_INPUT }
                PON_VEP = VCFANNO.out.vcf_enigma
                ch_versions = ch_versions.mix(VCFANNO.out.versions)
            }
            else {
                PON_VEP = VCFANNO_INPUT
            }

        // NON-OPTIONAL, NEEDED BY COYOTE
        ANNOTATE_VEP { PON_VEP }
        ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)

        MARK_GERMLINES { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(MARK_GERMLINES.out.versions)

        // BAF for CNVkit //
        FILTER_FOR_CNV { ANNOTATE_VEP.out.vcf_vep.join(concat_vcfs.filter { item -> item[1] == 'freebayes' })  }
        ch_versions = ch_versions.mix(FILTER_FOR_CNV.out.versions)

    emit:
        germline_variants   =   FILTER_FOR_CNV.out.vcf_only_germline    // channel: [ val(group), val(vc), file(vcf.gz) ]
        finished_vcf        =   MARK_GERMLINES.out.vcf_germline         // channel: [ val(group), val(vc), file(vcf.gz) ]
        versions            =   ch_versions                             // channel: [ file(versions) ]

}