#!/usr/bin/env nextflow

include { PON_FILTER                                } from '../../modules/local/filters/main'
include { FFPE_PON_FILTER                           } from '../../modules/local/filters/main'
include { ANNOTATE_VEP                              } from '../../modules/local/vep/main'
include { FIX_VEP_GNOMAD                            } from '../../modules/local/filters/main'
include { MARK_GERMLINES                            } from '../../modules/local/filters/main'
include { GERMLINE_FOR_CNVKIT                       } from '../../modules/local/filters/main'
include { VCFANNO                                   } from '../../modules/local/vcfanno/main'
include { POST_ANNOTATION_FILTERS                   } from '../../modules/local/filters/main'
include { BEDTOOLS_INTERSECT as CNV_BACKBONE_FILTER } from '../../modules/local/bedtools/main'
include { BEDTOOLS_INTERSECT as FILTER_FOR_CNV      } from '../../modules/local/bedtools/main'
include { TABIX_BGZIPTABIX as FILTER_FOR_CNV_TABIX  } from '../../modules/local/tabix/main'

workflow SNV_ANNOTATE {
    take: 
        agg_vcf         // channel: [mandatory] [ val(group), val(meta), file(agg.vcf) ]
        concat_vcfs     // channel: [mandatory] [ val(group), val(meta), val(vc), file(vcf.gz), file(vcf.gz.tbi) ]

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

        FIX_VEP_GNOMAD { ANNOTATE_VEP.out.vcf_vep }
        ch_versions = ch_versions.mix(FIX_VEP_GNOMAD.out.versions)

        MARK_GERMLINES { FIX_VEP_GNOMAD.out.fixed_vcf }
        ch_versions = ch_versions.mix(MARK_GERMLINES.out.versions)

        if (params.panel_profile_name.equals('gmshem') || params.panel_profile_name.equals('solid')) {
            if (!params.regions_backbone_idsnps) {
                error "params.regions_backbone_idsnps must be set for ${params.panel_profile_name} SNV annotation"
            }

            // CNV BACKBONE FILTER
            ch_cnv_backbone_input = MARK_GERMLINES.out.vcf_germline.map { group, meta, vcf ->
                tuple(group, meta, 'cnv_backbone', vcf)
            }
            ch_regions_backbone_idsnps = Channel.value(file(params.regions_backbone_idsnps))

            CNV_BACKBONE_FILTER ( ch_cnv_backbone_input, ch_regions_backbone_idsnps )
            ch_versions = ch_versions.mix( CNV_BACKBONE_FILTER.out.versions )
            POST_ANNOTATION_FILTERS ( CNV_BACKBONE_FILTER.out.intersected )
        } else {
            POST_ANNOTATION_FILTERS ( MARK_GERMLINES.out.vcf_germline )
        }

        ch_versions = ch_versions.mix(POST_ANNOTATION_FILTERS.out.versions)

        // BAF for CNVkit //
        ch_freebayes_vcfs = concat_vcfs
            .filter { group, meta, vc, vcf, tbi -> vc == 'freebayes' }
            .map { group, meta, vc, vcf, tbi -> [ group, vc, vcf ] }

        GERMLINE_FOR_CNVKIT ( ANNOTATE_VEP.out.vcf_vep )
        ch_versions = ch_versions.mix(GERMLINE_FOR_CNVKIT.out.versions)

        ch_filter_for_cnv_joined = ch_freebayes_vcfs.join(GERMLINE_FOR_CNVKIT.out.germline_vcf)

        ch_filter_for_cnv = ch_filter_for_cnv_joined
            .map { group, vc, vcf_unfilt, meta, germline_vcf ->
                [ group, meta, vc, vcf_unfilt ]
            }

        ch_filter_for_cnv_germline = ch_filter_for_cnv_joined.map { group, vc, vcf_unfilt, meta, germline_vcf -> germline_vcf }

        FILTER_FOR_CNV ( ch_filter_for_cnv, ch_filter_for_cnv_germline )
        ch_versions = ch_versions.mix(FILTER_FOR_CNV.out.versions)

        FILTER_FOR_CNV_TABIX ( FILTER_FOR_CNV.out.vcf_intersected )
        ch_versions = ch_versions.mix(FILTER_FOR_CNV_TABIX.out.versions)

    emit:
        germline_variants   =   FILTER_FOR_CNV_TABIX.out.gz_index.map { group, meta, vc, vcf, tbi -> [ group, meta, vcf, tbi ] } // channel: [ val(group), val(meta), file(vcf.gz), file(vcf.gz.tbi) ]
        finished_vcf        =   POST_ANNOTATION_FILTERS.out.filtered_vcf         // channel: [ val(group), val(vc), file(vcf.gz) ]
        vep_vcf             =   ANNOTATE_VEP.out.vcf_vep                         // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        versions            =   ch_versions                                      // channel: [ file(versions) ]
}
