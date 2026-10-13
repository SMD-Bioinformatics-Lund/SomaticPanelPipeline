#!/usr/bin/env nextflow

include { ANNOTATE_VEP                          				} from '../../modules/local/vep/main'
include { COYOTE_SEGMENTS                       				} from '../../modules/local/filters/main'
include { MERGE_SEGMENTS                        				} from '../../modules/local/filters/main'
include { BEDTOOLS_INTERSECT as GENE_INTERSECT  				} from '../../modules/local/bedtools/main'
include { COYOTE_SEGMENTS_JSON                  				} from '../../modules/local/filters/main'
include { MERGE_JSON                            				} from '../../modules/local/filters/main'
include { SVDB_QUERY as SVDB_QUERY_ALL                        } from '../../modules/local/svdb/main'
include { SVDB_QUERY as SVDB_QUERY_GC                         } from '../../modules/local/svdb/main'
include { SVDB_QUERY as SVDB_QUERY_AT                         } from '../../modules/local/svdb/main'

workflow CNV_ANNOTATE {
	take: 
		tumor              // channel: [mandatory] [ val(group), val(meta), file(tumor_merged_vcf) ]
		normal             // channel: [mandatory] [ val(group), val(meta), file(normal_merged_vcf) ]
		// meta               // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]
	main:
		ch_versions = Channel.empty()

		if (params.loqusdb_export) {
			def databases = params.loqusdb_vcfs
			if (!(databases instanceof Map) || databases.isEmpty()) {
				error 'SVDB annotation requires a non-empty database map using ALL, GC and/or AT'
			}
			if (databases.keySet().toList() != ['ALL', 'GC', 'AT'].findAll { databases.containsKey(it) }) {
				error 'SVDB databases must use ALL, GC and/or AT keys in that order'
			}
			annotated_vcf = tumor
			if (databases.containsKey('ALL')) {
				SVDB_QUERY_ALL(annotated_vcf, Channel.value(tuple([id: 'ALL'], file(databases.ALL))))
				ch_versions = ch_versions.mix(SVDB_QUERY_ALL.out.versions)
				annotated_vcf = SVDB_QUERY_ALL.out.queried_vcf
			}
			if (databases.containsKey('GC')) {
				SVDB_QUERY_GC(annotated_vcf, Channel.value(tuple([id: 'GC'], file(databases.GC))))
				ch_versions = ch_versions.mix(SVDB_QUERY_GC.out.versions)
				annotated_vcf = SVDB_QUERY_GC.out.queried_vcf
			}
			if (databases.containsKey('AT')) {
				SVDB_QUERY_AT(annotated_vcf, Channel.value(tuple([id: 'AT'], file(databases.AT))))
				ch_versions = ch_versions.mix(SVDB_QUERY_AT.out.versions)
				annotated_vcf = SVDB_QUERY_AT.out.queried_vcf
			}
			gene_input = annotated_vcf.mix(normal)
		}

		else {
			gene_input = tumor.mix(normal)
		}

		ch_gene_intersect_input = gene_input.map { group, meta, vcf ->
			tuple(group, meta, 'genes', vcf)
		}
		ch_gene_gtf = Channel.value(file(params.gene_gtf))

		GENE_INTERSECT ( ch_gene_intersect_input.combine(ch_gene_gtf) )
		ch_versions = ch_versions.mix(GENE_INTERSECT.out.versions)
		COYOTE_SEGMENTS_JSON ( GENE_INTERSECT.out.intersected )
		ch_versions = ch_versions.mix(COYOTE_SEGMENTS_JSON.out.versions)

		COYOTE_SEGMENTS ( tumor.mix(normal) )

		MERGE_SEGMENTS ( COYOTE_SEGMENTS.out.filtered.groupTuple() )
		ch_versions = ch_versions.mix(MERGE_SEGMENTS.out.versions)
		MERGE_JSON ( COYOTE_SEGMENTS_JSON.out.json_panel.groupTuple() )
		ch_versions = ch_versions.mix(MERGE_JSON.out.versions)
		COYOTE_JSON = MERGE_JSON.out.merged

		ch_versions = ch_versions.mix(COYOTE_SEGMENTS.out.versions)

	emit:
		segments 	= 	MERGE_SEGMENTS.out.merged	// channel: [ val(group), file(cn-segments.panel.merged.bed) ]
		s_json      =   MERGE_JSON.out.merged       // channel: [ val(group), file(panel.json) ]
		versions    =   ch_versions 				// channel: [ file(versions) ]
}
