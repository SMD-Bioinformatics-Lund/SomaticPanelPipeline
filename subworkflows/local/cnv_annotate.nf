#!/usr/bin/env nextflow

include { ANNOTATE_VEP                          				} from '../../modules/local/vep/main'
include { COYOTE_SEGMENTS_BED                   				} from '../../modules/local/filters/main'
include { COYOTE_SEGMENTS                       				} from '../../modules/local/filters/main'
include { MERGE_SEGMENTS                        				} from '../../modules/local/filters/main'
include { BEDTOOLS_INTERSECT as GENE_INTERSECT  				} from '../../modules/local/bedtools/main'
include { BEDTOOLS_INTERSECT as COYOTE_SEGMENTS_INTERSECT 		} from '../../modules/local/bedtools/main'
include { COYOTE_SEGMENTS_JSON                  				} from '../../modules/local/filters/main'
include { MERGE_JSON                            				} from '../../modules/local/filters/main'
include { SVDB_ANNOTATE_ARTEFACTS               				} from '../../modules/local/svdb/main'

workflow CNV_ANNOTATE {
	take: 
		tumor              // channel: [mandatory] [ val(group), val(meta), file(tumor_merged_vcf) ]
		normal             // channel: [mandatory] [ val(group), val(meta), file(normal_merged_vcf) ]
		// meta               // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
	main:
		ch_versions = Channel.empty()

		if (params.loqusdb_export) {
			SVDB_ANNOTATE_ARTEFACTS(tumor, params.loqusdb_vcfs)
			gene_input = SVDB_ANNOTATE_ARTEFACTS.out.artefacts.mix(normal)
		}

		else {
			gene_input = tumor.mix(normal)
		}

		ch_gene_intersect_input = gene_input.map { group, meta, vcf ->
			tuple(group, meta, 'genes', vcf)
		}
		ch_gene_gtf = Channel.value(file(params.gene_gtf))

		GENE_INTERSECT ( ch_gene_intersect_input, ch_gene_gtf )
		COYOTE_SEGMENTS_JSON ( GENE_INTERSECT.out.intersected )

		COYOTE_SEGMENTS_BED ( tumor.mix(normal) )
		ch_gencode_genes = Channel.value(file(params.gencode_genes))
		COYOTE_SEGMENTS_INTERSECT ( COYOTE_SEGMENTS_BED.out.raw_segments_bed, ch_gencode_genes )
		COYOTE_SEGMENTS ( COYOTE_SEGMENTS_INTERSECT.out.intersected )

		MERGE_SEGMENTS ( COYOTE_SEGMENTS.out.filtered.groupTuple() )
		MERGE_JSON ( COYOTE_SEGMENTS_JSON.out.json_panel.groupTuple() )
		COYOTE_JSON = MERGE_JSON.out.merged

		ch_versions = ch_versions.mix(COYOTE_SEGMENTS_BED.out.versions)
		ch_versions = ch_versions.mix(COYOTE_SEGMENTS_INTERSECT.out.versions)
		ch_versions = ch_versions.mix(COYOTE_SEGMENTS.out.versions)

	emit:
		segments 	= 	MERGE_SEGMENTS.out.merged	// channel: [ val(group), file(cn-segments.panel.merged.bed) ]
		s_json      =   MERGE_JSON.out.merged       // channel: [ val(group), file(panel.json) ]
		versions    =   ch_versions 				// channel: [ file(versions) ]
}
