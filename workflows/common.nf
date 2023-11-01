#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { SNV_ANNOTATE                  } from '../subworkflows/local/snv_annotate'
include { CNV_CALLING                   } from '../subworkflows/local/cnv_calling'
include { BIOMARKERS                    } from '../subworkflows/local/biomarkers'
include { QC                            } from '../subworkflows/local/qc'
include { ADD_TO_DB                     } from '../subworkflows/local/add_to_db'
include { CNV_ANNOTATE                  } from '../subworkflows/local/cnv_annotate'
include { FUSIONS                       } from '../subworkflows/local/fusions'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'

csv = file(params.csv)

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


workflow SPP_COMMON {

    ch_versions = Channel.empty()

    // Checks input, creates meta-channel and decides whether data should be downsampled //
    CHECK_INPUT ( Channel.fromPath(csv) )

    // Downsample if meta.sub == value and not false //
    SAMPLE ( CHECK_INPUT.out.fastq )  
    .set{ ch_trim }
    ch_versions = ch_versions.mix(ch_trim.versions)

    // Do alignment if downsample was false and mix with SAMPLE subworkflow output
    ALIGN_SENTIEON ( 
        ch_trim.fastq_trim,
        CHECK_INPUT.out.meta
    )
    .set { ch_mapped }
    ch_versions = ch_versions.mix(ch_mapped.versions)

    QC ( 
        ch_mapped.qc_out, 
        ch_mapped.bam_lowcov
    )
    .set { ch_qc }
    ch_versions = ch_versions.mix(ch_qc.versions)

    SNV_CALLING ( 
        ch_mapped.bam_umi.groupTuple(),
        ch_mapped.bam_dedup,
        beds,
        CHECK_INPUT.out.meta,
        ch_qc.melt_qc
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

    CNV_CALLING ( 
        ch_mapped.bam_umi, 
        ch_vcf_anno.germline_variants,
        CHECK_INPUT.out.meta,
        ch_mapped.bam_dedup,
        gatk_ref
    )
    .set { ch_cnvcalled }
    ch_versions = ch_versions.mix(ch_cnvcalled.versions)

    CNV_ANNOTATE (
        ch_cnvcalled.tumor_vcf,
        ch_cnvcalled.normal_vcf,
        CHECK_INPUT.out.meta
    )
    .set { ch_cnv }
    ch_versions = ch_versions.mix(ch_cnv.versions)


    FUSIONS (
        ch_trim.fastq_trim,
        CHECK_INPUT.out.meta,
        ch_mapped.bam_dedup
    )
    .set { ch_fusion }
    ch_versions = ch_versions.mix(ch_fusion.versions)

    BIOMARKERS (
        CHECK_INPUT.out.meta,
        ch_cnvcalled.cnvkit_hrd,
        ch_mapped.bam_umi, 
        ch_mapped.bam_dedup
    )
    .set { ch_bio }
    ch_versions = ch_versions.mix(ch_bio.versions)


    ADD_TO_DB (
        ch_vcf_anno.finished_vcf,
        ch_qc.lowcov.filter { item -> item[1] == 'T' },
        ch_cnv.segments,
        ch_cnvcalled.gens,
        ch_cnvcalled.gatcov_plot,
        ch_fusion.fusions,
        ch_bio.biomarkers
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
