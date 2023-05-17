#!/usr/bin/env nextflow

include { MANTA as MANTA_FUSIONS               } from '../../modules/local/manta/main'
include { SVDB_MERGE_PANEL as JOIN_FUSIONS     } from '../../modules/local/svdb/main'
include { GENEFUSE                             } from '../../modules/local/genefuse/main'
include { SNPEFF                               } from '../../modules/local/snpeff/main'
include { FILTER_MANTA                         } from '../../modules/local/filters/main'
include { GENEFUSE_JSON_TO_VCF                 } from '../../modules/local/filters/main'

workflow FUSIONS {
	take: 
        fastq                // val(group), val(meta), file(r1), file(r2)
		meta                 // map: (csv meta info)
		bam_markdup          // val(group), val(meta), file(bam), file(bai), file(bqsr) markdup bam

    main:
        // genefuse //
        GENEFUSE(fastq.filter { it -> it[1].type == "T" })
        GENEFUSE_JSON_TO_VCF(GENEFUSE.out.genefuse_json)
        // join meta-info. fastq-meta-channel differs from global meta-channel
        GENEFUSE_TUMOR = GENEFUSE_JSON_TO_VCF.out.genefuse_vcf.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }
        // manta //
        MANTA_FUSIONS(bam_markdup.filter { it -> it[1].type == "T" }.groupTuple(), params.mantafusions, "fusions")
        MANTA_FUSION_TUMOR = MANTA_FUSIONS.out.manta_vcf_tumor.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }
        FILTER_MANTA(MANTA_FUSION_TUMOR)
        JOIN_FUSIONS( GENEFUSE_TUMOR.mix(FILTER_MANTA.out.bnd_filtered).groupTuple(by:[0,1]) )
        SNPEFF(JOIN_FUSIONS.out.merged_vcf)

        

    emit:
        fusions = SNPEFF.out.snpeff_vcf

}