#!/usr/bin/env nextflow

include { MANTA as MANTA_FUSIONS               } from '../../modules/local/manta/main'
include { SVDB_MERGE_PANEL as JOIN_FUSIONS     } from '../../modules/local/svdb/main'
include { GENEFUSE                             } from '../../modules/local/genefuse/main'
include { SNPEFF                               } from '../../modules/local/snpeff/main'
include { FILTER_MANTA                         } from '../../modules/local/filters/main'
include { GENEFUSE_JSON_TO_VCF                 } from '../../modules/local/filters/main'

workflow FUSIONS {
	take: 
        fastq                // channel: [mandatory] [ val(group), val(meta), file(r1), file(r2) ]
		// meta                 // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
		bam_markdup          // channel: [mandatory] [ val(group), val(meta), file(marked_bam), file(bai), file(bqsr), file(dedup_metrics.txt) ]

    main:
        ch_versions = Channel.empty()

        // genefuse //
        GENEFUSE(fastq.filter { it -> it[1].type == "T" })
        ch_versions = ch_versions.mix(GENEFUSE.out.versions)

        GENEFUSE_JSON_TO_VCF(GENEFUSE.out.genefuse_json)
        ch_versions = ch_versions.mix(GENEFUSE_JSON_TO_VCF.out.versions)

        // manta //
        MANTA_FUSIONS(bam_markdup.filter { it -> it[1].type == "T" }.groupTuple(), params.mantafusions, "fusions")
        ch_versions         = ch_versions.mix(MANTA_FUSIONS.out.versions)
        MANTA_FUSION_TUMOR  = MANTA_FUSIONS.out.manta_vcf.flatMap { group, meta, vcfs ->
            files = vcfs instanceof List ? vcfs : [vcfs]
            tumor_meta = meta.find { it.type == "T" || it.type == "tumor" }
            tumor_meta ? files.findAll { it.name.startsWith("${tumor_meta.id}_manta.") }.collect { vcf -> tuple(group, tumor_meta, vcf) } : []
        }

        FILTER_MANTA(MANTA_FUSION_TUMOR)
        ch_versions = ch_versions.mix(FILTER_MANTA.out.versions)

        GENEFUSE_TUMOR_LABELLED = GENEFUSE_JSON_TO_VCF.out.genefuse_vcf.map { group, meta, vcf ->
            tuple(group, meta, 'genefuse', vcf)
        }
        MANTA_FUSION_TUMOR_LABELLED = FILTER_MANTA.out.bnd_filtered.map { group, meta, vcf ->
            tuple(group, meta, 'manta', vcf)
        }

        JOIN_FUSIONS( GENEFUSE_TUMOR_LABELLED.mix(MANTA_FUSION_TUMOR_LABELLED).groupTuple(by:[0,1]) )
        ch_versions = ch_versions.mix(JOIN_FUSIONS.out.versions)

        SNPEFF(JOIN_FUSIONS.out.merged_vcf)
        ch_versions = ch_versions.mix(SNPEFF.out.versions)
        output = Channel.empty().mix(SNPEFF.out.snpeff_vcf)

    emit:
        fusions     =   output          // channel: [ val(group), file(merged.annotated.vcf) ]
        versions    =   ch_versions     // channel: [ file(versions) ]

}
