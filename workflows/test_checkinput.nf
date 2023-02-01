#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { CNV_CALLING                   } from '../subworkflows/local/cnv_calling'
include { BIOMARKERS                    } from '../subworkflows/local/biomarkers'
include { QC                            } from '../subworkflows/local/qc'
include { ADD_TO_DB                     } from '../subworkflows/local/add_to_db'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { CNV_ANNOTATE                  } from '../subworkflows/local/cnv_annotate'

println(params.genome_file)

csv = file(params.csv)
println(csv)

cnvkit_reference_exons   = params.cnvkit_reference_exons      ? Channel.fromPath(params.cnvkit_reference_exons).map {it -> [[id:it[0].simpleName], it]}.collect(): ( ch_references.cnvkit_reference_exons                ?: Channel.empty() )
cnvkit_reference_backbone   = params.cnvkit_reference_backbone      ? Channel.fromPath(params.cnvkit_reference_backbone).map {it -> [[id:it[0].simpleName], it]}.collect(): ( ch_references.cnvkit_reference_backbone                ?: Channel.empty() )

// Split bed file in to smaller parts to be used for parallel variant calling
Channel
    .fromPath("${params.regions_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.regions_bed}" }
    .splitText( by: 1000, file: 'bedpart.bed' )
    .set { beds }

Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.set{ gatk_ref}


workflow SOLID_GMS {

	// Checks input, creates meta-channel and decides whether data should be downsampled //
	CHECK_INPUT ( Channel.fromPath(csv) )

	// Downsample if meta.sub == value and not false //
	SAMPLE ( CHECK_INPUT.out.fastq )  
	.set{ ch_trim }
	// Do alignment if downsample was false and mix with SAMPLE subworkflow output
	ALIGN_SENTIEON ( 
		ch_trim.fastq_trim,
		CHECK_INPUT.out.meta
	)
	.set { ch_mapped } 
	QC ( 
        ch_mapped.qc_out, 
		ch_mapped.bam_lowcov
    )
	.set { ch_qc }
	SNV_CALLING ( 
		ch_mapped.bam_umi.groupTuple(),
		beds,
		CHECK_INPUT.out.meta
	)
	.set { ch_vcf }
	CNV_CALLING ( 
		ch_mapped.bam_umi, 
		ch_vcf.germline_variants,
		CHECK_INPUT.out.meta,
		ch_mapped.bam_dedup,
		gatk_ref
	)
	.set { ch_cnvcalled }
	CNV_ANNOTATE (
		ch_cnvcalled.tumor_vcf,
		ch_cnvcalled.normal_vcf,
		CHECK_INPUT.out.meta
	)
	.set { ch_cnv }
	ADD_TO_DB (
		ch_vcf.finished_vcf,
		ch_qc.lowcov.filter { item -> item[1] == 'T' },
		ch_cnv.segments,
		ch_cnvcalled.gens,
		ch_cnvcalled.gatcov_plot
	)


}

workflow {
	SOLID_GMS()
}













workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}