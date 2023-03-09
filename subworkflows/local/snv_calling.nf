#!/usr/bin/env nextflow

include { FREEBAYES                } from '../../modules/local/freebayes/main'
include { VARDICT                  } from '../../modules/local/vardict/main'
include { CONCATENATE_VCFS         } from '../../modules/local/concatenate_vcfs/main'
include { AGGREGATE_VCFS           } from '../../modules/local/concatenate_vcfs/main'
include { PON_FILTER               } from '../../modules/local/filters/main'
include { FFPE_PON_FILTER          } from '../../modules/local/filters/main'
include { ANNOTATE_VEP             } from '../../modules/local/filters/main'
include { MARK_GERMLINES           } from '../../modules/local/filters/main'
include { FILTER_FOR_CNV           } from '../../modules/local/filters/main'
include { CONTAMINATION            } from '../../modules/local/concatenate_vcfs/main'

workflow SNV_CALLING {
    take: 
        bam_umi
        beds
        meta

    main:
        // Variantcallers //
        // split by bed-file to speed up calling //
        FREEBAYES ( bam_umi, beds)
        VARDICT ( bam_umi, beds)
        // Prepare vcf parts for concatenation //
        vcfparts_freebayes = FREEBAYES.out.vcfparts_freebayes.groupTuple(by:[0,1])
        vcfparts_vardict   = VARDICT.out.vcfparts_vardict.groupTuple(by:[0,1])
        //vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)
        vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict)
        // Join vcfs split by bedparts //
        CONCATENATE_VCFS { vcfs_to_concat }
        // Aggregate all callers to one VCF
        AGGREGATE_VCFS { CONCATENATE_VCFS.out.concatenated_vcfs.groupTuple().join(meta.groupTuple()) }
        // Filter with PoN, annotate with VEP, mark germlines
        PON_FILTER { AGGREGATE_VCFS.out.vcf_concat }
        if (params.assay == "solid") {
            FFPE_PON_FILTER { PON_FILTER.out.vcf_pon } // needs to be merged with normal PON above, myeloid should not be be FFPE annotated
        }
        
        ANNOTATE_VEP { FFPE_PON_FILTER.out.vcf_pon_ffpe } 
        MARK_GERMLINES { ANNOTATE_VEP.out.vcf_vep }
        // filter for CNVkit //
        FILTER_FOR_CNV { ANNOTATE_VEP.out.vcf_vep.join(CONCATENATE_VCFS.out.concatenated_vcfs.filter { item -> item[1] == 'freebayes' })  }
        // contamination values from VCF //
        CONTAMINATION { ANNOTATE_VEP.out.vcf_vep }

    emit:
        concat_vcfs = CONCATENATE_VCFS.out.concatenated_vcfs
        germline_variants = FILTER_FOR_CNV.out.vcf_only_germline
        finished_vcf = MARK_GERMLINES.out.vcf_germline

}