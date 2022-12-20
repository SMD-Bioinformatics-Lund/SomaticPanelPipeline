#!/usr/bin/env nextflow

include { CNVKIT                               } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH                         } from '../../modules/local/cnvkit/main'
include { CNVKIT_GENS                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_PLOT                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_CALL                          } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_BACKBONE      } from '../../modules/local/cnvkit/main'
include { CNVKIT_BATCH as CNVKIT_EXONS         } from '../../modules/local/cnvkit/main'
include { GATKCOV_BAF                          } from '../../modules/local/GATK/main'
include { GATKCOV_COUNT                        } from '../../modules/local/GATK/main'
include { GATKCOV_CALL                         } from '../../modules/local/GATK/main'


workflow CNV_CALLING {
	take: 
		bam_umi              // val(group), val(meta), file(bam), file(bai), file(bqsr)
		germline_variants    // val(group), file(vcf), file(tbi)
		meta                 // map: (csv meta info)

	main:
		////////////////////////// CNVKIT ////////////////////////////////////////////////////
		// if backbone + exon pool differs in pool ratio do backbone and exons separatly
		if (!params.cnvkit_split) {
			CNVKIT_BATCH ( bam_umi, params.cnvkit_reference )
			batch_plot_cns = CNVKIT_BATCH.out.cnvkit_cns
			batch_plot_cnr = CNVKIT_BATCH.out.cnvkit_cnr
		}
		else {
			CNVKIT_EXONS ( bam_umi, params.cnvkit_reference_exons )
			batch_plot_cnr = CNVKIT_EXONS.out.cnvkit_cnr
			batch_plot_cns = CNVKIT_EXONS.out.cnvkit_cns
			CNVKIT_BACKBONE ( bam_umi, params.cnvkit_reference_backbone )
		}
		// call, plot and export segments ::: cnvkit
		CNVKIT_PLOT ( batch_plot_cns.join(batch_plot_cnr, by:[0,1]).join(germline_variants) )
		CNVKIT_CALL ( batch_plot_cns.join(batch_plot_cnr, by:[0,1]).join(germline_variants) )
		CNVKIT_GENS ( batch_plot_cnr.join(germline_variants) )
		///////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////// GATK SEGMENT CALLING /////////////////////////////////////
		GATKCOV_BAF ( bam_umi )
		GATKCOV_COUNT ( bam_umi )
		GATKCOV_CALL { GATKCOV_BAF.out.gatk_baf.join(GATKCOV_COUNT.out.gatk_count, by:[0,1]) }
		///////////////////////////////////////////////////////////////////////////////////////
		
//CNVKIT { bam_umi.join(meta, by:[0,1,2]).join(germline_variants) }
//CNVKIT { bam_umi.join(meta, by:[0,1,2]).join(concat_vcf.filter { item -> item[1] == 'freebayes' } ).view() }
	emit:
		// baflogr = CNVKIT_.out.baflogr
		// cnvkit_cns = CNVKIT_BATCH.out.cnvkit_cns
		// cnvkitsegment_purity = CNVKIT.out.cnvkitsegment_purity
		test = GATKCOV_BAF.out.gatk_baf

}