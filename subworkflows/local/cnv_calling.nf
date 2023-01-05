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
include { GATK2VCF                             } from '../../modules/local/GATK/main'
include { GATK_COUNT_GERMLINE                  } from '../../modules/local/GATK/main'
include { GATK_CALL_PLOIDY                     } from '../../modules/local/GATK/main'
include { GATK_CALL_GERMLINE_CNV               } from '../../modules/local/GATK/main'
include { FILTER_MERGE_GATK                    } from '../../modules/local/GATK/main'
include { POSTPROCESS                          } from '../../modules/local/GATK/main'
include { MANTA                                } from '../../modules/local/manta/main'
include { SVDB_MERGE_PANEL as JOIN_TUMOR       } from '../../modules/local/svdb/main'
include { SVDB_MERGE_PANEL as JOIN_NORMAL      } from '../../modules/local/svdb/main'


workflow CNV_CALLING {
	take: 
		bam_umi              // val(group), val(meta), file(bam), file(bai), file(bqsr) umi deduped bam
		germline_variants    // val(group), file(vcf), file(tbi)
		meta                 // map: (csv meta info)
		bam_markdup          // val(group), val(meta), file(bam), file(bai), file(bqsr) markdup bam
		gatk_ref             // val(interger), val(part_of_genome) used for germline gatk-calling

	main:
		////////////////////////// CNVKIT /////////////////////////////////////////////////////
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
			cnvkit_hrd = CNVKIT_CALL.out.cnvkitsegment
		}
		else {
			CNVKIT_BATCH ( bam_umi, params.cnvkit_reference, "full" )
			CNVKIT_EXONS ( bam_umi, params.cnvkit_reference_exons, "exons" )
			CNVKIT_BACKBONE ( bam_umi, params.cnvkit_reference_backbone, "backbone" )
			// call, plot and export segments ::: cnvkit
			CNVKIT_PLOT ( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).combine(germline_variants, by:[0]) )
			//CNVKIT_CALL ( CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3]).mix(CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3]),CNVKIT_BATCH.out.cnvkit_cns.join(CNVKIT_BATCH.out.cnvkit_cnr, by:[0,1,3])).combine(germline_variants, by:[0]) )
			CNVKIT_CALL ( CNVKIT_EXONS.out.cnvkit_cns.join(CNVKIT_EXONS.out.cnvkit_cnr, by:[0,1,3]).mix(CNVKIT_BACKBONE.out.cnvkit_cns.join(CNVKIT_BACKBONE.out.cnvkit_cnr, by:[0,1,3])).combine(germline_variants, by:[0]) )
			CNVKIT_GENS ( CNVKIT_EXONS.out.cnvkit_cnr.mix(CNVKIT_BACKBONE.out.cnvkit_cnr).combine(germline_variants, by:[0]) )
			MERGE_GENS  ( CNVKIT_GENS.out.cnvkit_gens.groupTuple(by:[0,1]) )
			cnvkitplot = CNVKIT_PLOT.out.cnvkitplot.filter { it -> it[2] == "backbone" }
			cnvkit_hrd = CNVKIT_CALL.out.cnvkitsegment.filter { it -> it[2] == "backbone" }
		}

		///////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
		// Do calling somatic CNV calling, use normal allelic counts for somatic as well     //
		///////////////////////////////////////////////////////////////////////////////////////
		GATKCOV_BAF ( bam_umi )
		GATKCOV_COUNT ( bam_umi )
		GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count,by:[0,1]).groupTuple() }
		GATK2VCF ( GATKCOV_CALL.out.gatcov_called.join(meta.filter( it -> it[1].type == "T" )) )
		// Do germline calling for normal
		GATK_COUNT_GERMLINE ( bam_umi.filter { it -> it[1].type == "N" })
		GATK_CALL_PLOIDY ( GATK_COUNT_GERMLINE.out.count_germline )
		GATK_CALL_GERMLINE_CNV( GATK_COUNT_GERMLINE.out.count_germline.join(GATK_CALL_PLOIDY.out.gatk_ploidy,by:[0,1]).groupTuple(by:[0,1]).combine(gatk_ref) )
		CALLED = GATK_CALL_GERMLINE_CNV.out.gatk_call_germline.groupTuple(by:[0,1])
		PLOIDY = GATK_CALL_PLOIDY.out.gatk_ploidy.groupTuple(by:[0,1])
		POSTPROCESS ( CALLED.join(PLOIDY,by:[0,1]).combine(gatk_ref.groupTuple(by:[3])))
		FILTER_MERGE_GATK ( POSTPROCESS.out.gatk_germline_segmentsvcf )
		/////////////////////////// MANTA /////////////////////////////////////////////////////
		MANTA ( bam_markdup.groupTuple() )

		// Join germline vcf
		MANTA_NORMAL = MANTA.out.manta_vcf_normal.join(meta.filter( it -> it[1].type == "N" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }
		GATK_NORMAL = FILTER_MERGE_GATK.out.gatk_normal_vcf.join(meta.filter( it -> it[1].type == "N" ) ).map{ val-> tuple(val[0], val[2], val[1] )}
		JOIN_NORMAL ( GATK_NORMAL.mix(MANTA_NORMAL).groupTuple(by:[0,1]) )
		// Join tumor vcf
		GATK_TUMOR = GATK2VCF.out.tumor_vcf
		MANTA_TUMOR = MANTA.out.manta_vcf_tumor.join(meta.filter( it -> it[1].type == "T" ) ).map{ val-> tuple(val[0], val[2], val[1] ) }
		JOIN_TUMOR ( GATK_TUMOR.mix(MANTA_TUMOR).groupTuple(by:[0,1]) )
		

	emit:
		cnvkit_plot = cnvkitplot
		cnvkit_hrd = cnvkit_hrd
		tumor_vcf = JOIN_TUMOR.out.merged_vcf
		normal_vcf = JOIN_NORMAL.out.merged_vcf

}