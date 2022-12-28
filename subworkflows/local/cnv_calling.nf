#!/usr/bin/env nextflow

include { CNVKIT                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH                         } from '../../modules/local/cnvkit/main'
include { CNVKIT_GENS                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_PLOT                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL                          } from '../../modules/local/cnvkit/main'
include { MERGE_GENS                           } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_BACKBONE      } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_EXONS         } from '../../modules/local/cnvkit/main'
include { GATKCOV_BAF                          } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT                        } from '../../modules/local/GATK/main'
include { GATKCOV_CALL                         } from '../../modules/local/GATK/main'
include { MANTA                                } from '../../modules/local/manta/main'


workflow CNV_CALLING {
	take: 
		bam_umi              // val(group), val(meta), file(bam), file(bai), file(bqsr) umi deduped bam
		germline_variants    // val(group), file(vcf), file(tbi)
		meta                 // map: (csv meta info)
		bam_markdup          // val(group), val(meta), file(bam), file(bai), file(bqsr) markdup bam

	main:
		////////////////////////// CNVKIT ////////////////////////////////////////////////////
		// if backbone + exon pool differs in pool ratio do backbone and exons separatly
		if (!params.cnvkit_split) {
			CNVKIT_BATCH ( bam_umi, params.cnvkit_reference, "full" )
			batch_plot_cns = CNVKIT_BATCH.out.cnvkit_cns
			batch_plot_cnr = CNVKIT_BATCH.out.cnvkit_cnr
			// call, plot and export segments ::: cnvkit
			CNVKIT_PLOT ( batch_plot_cns.join(batch_plot_cnr, by:[0,1,3]).combine(germline_variants, by:[0]) )
			CNVKIT_CALL ( batch_plot_cns.join(batch_plot_cnr, by:[0,1,3]).combine(germline_variants, by:[0]) )
			CNVKIT_GENS ( batch_plot_cnr.mix(CNVKIT_BACKBONE.out.cnvkit_cnr).combine(germline_variants, by:[0]) )
			MERGE_GENS  ( CNVKIT_GENS.out.cnvkit_gens.groupTuple(by:[0,1]) )
			cnvkitplot = CNVKIT_PLOT.out.cnvkitplot
		}
		else {
			CNVKIT_BATCH ( bam_umi, params.cnvkit_reference, "full" )
			CNVKIT_EXONS ( bam_umi, params.cnvkit_reference_exons, "exons" )
			CNVKIT_BACKBONE ( bam_umi, params.cnvkit_reference_backbone, "backbone" )
			// call, plot and export segments ::: cnvkit
			CNVKIT_PLOT ( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).mix(CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3]),CNVKIT_BATCH.out.cnvkit_cns.join(CNVKIT_BATCH.out.cnvkit_cnr, by:[0,1,3])).combine(germline_variants, by:[0]) )
			CNVKIT_CALL ( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).mix(CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3]),CNVKIT_BATCH.out.cnvkit_cns.join(CNVKIT_BATCH.out.cnvkit_cnr, by:[0,1,3])).combine(germline_variants, by:[0]) )
			CNVKIT_GENS ( CNVKIT_EXONS.out.cnvkit_cnr.mix(CNVKIT_BACKBONE.out.cnvkit_cnr).combine(germline_variants, by:[0]) )
			MERGE_GENS  ( CNVKIT_GENS.out.cnvkit_gens.groupTuple(by:[0,1]) )
			cnvkitplot = CNVKIT_PLOT.out.cnvkitplot.filter { it -> it[2] == "backbone" }
		}

		///////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
		GATKCOV_BAF ( bam_umi )
		GATKCOV_COUNT ( bam_umi )
		GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count, by:[0,1]) } // maybe add normal allelic if paired.
		///////////////////////////////////////////////////////////////////////////////////////

		/////////////////////////// MANTA /////////////////////////////////////////////////////
		MANTA ( bam_markdup.groupTuple().view() )
		

	emit:
		cnvkit_plot = cnvkitplot
		test = GATKCOV_BAF.out.gatk_baf

}