#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { SNV_ANNOTATE                  } from '../subworkflows/local/snv_annotate'
include { BAM_QC                        } from '../subworkflows/local/bam_qc'
include { VCF_QC                        } from '../subworkflows/local/vcf_qc'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'

csv = file(params.csv)

params.paired = csv.countLines() > 2 ? true : false


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



workflow SPP_QC {

    ch_versions = Channel.empty()

    // Checks input, creates meta-channel and decides whether data should be downsampled //
    CHECK_INPUT ( Channel.fromPath(csv), params.paired )

    // Downsample if meta.sub == value and not false //
    SAMPLE ( CHECK_INPUT.out.fastq )  
    .set{ ch_trim }
    ch_versions = ch_versions.mix(ch_trim.versions)

    // Do alignment if downsample was false and mix with SAMPLE subworkflow output
    ALIGN_SENTIEON ( 
        ch_trim.fastq_trim,
        CHECK_INPUT.out.bam,
        CHECK_INPUT.out.meta
    )
    .set { ch_mapped }
    ch_versions = ch_versions.mix(ch_mapped.versions)

    BAM_QC (
        ch_mapped.bam_umi,
        ch_mapped.bam_dedup,
        ch_mapped.dedup_metrics
    )
    .set { ch_qc }
    ch_versions = ch_versions.mix(ch_qc.versions)

    SNV_CALLING ( 
        ch_mapped.bam_umi.groupTuple(),
        ch_mapped.bam_dedup,
        beds,
        CHECK_INPUT.out.meta,
        ch_qc.melt_qc,
        ch_qc.dedup_bam_is_metrics.groupTuple(),
    )
    .set { ch_vcf }
    ch_versions = ch_versions.mix(ch_vcf.versions)

    SNV_ANNOTATE (
        ch_vcf.agg_vcf,
        ch_vcf.concat_vcfs,
        CHECK_INPUT.out.meta
    )
    .set { ch_vcf_anno }
    ch_versions = ch_versions.mix(ch_vcf_anno.versions)

    VCF_QC (
        ch_vcf_anno.vep_vcf,
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        CHECK_INPUT.out.meta
    )

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
        csv         : ${params.csv}
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
    logFile = file("${params.resultsdir}/cron/logs/" + base + ".complete")
    logFile.text = msg
    logFile.append(error)
}