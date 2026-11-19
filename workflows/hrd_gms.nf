#!/usr/bin/env nextflow


nextflow.enable.dsl = 2


include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { BAM_QC                        } from '../subworkflows/local/bam_qc'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { SNV_ANNOTATE                  } from '../subworkflows/local/snv_annotate'
include { CNV_CALLING                   } from '../subworkflows/local/cnv_calling'
include { BIOMARKERS                    } from '../subworkflows/local/biomarkers'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'

println(params.genome_file)
genome_file = file(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

csv = file(params.csv)
println(csv)


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

    ch_versions = Channel.empty()

    CHECK_INPUT ( Channel.fromPath(csv) )

    ALIGN_SENTIEON ( CHECK_INPUT.out.fastq, CHECK_INPUT.out.bam, CHECK_INPUT.out.meta )
    .set { ch_mapped }
    ch_versions = ch_versions.mix(ch_mapped.versions)

    BAM_QC ( ch_mapped.bam_umi, ch_mapped.bam_dedup, ch_mapped.dedup_metrics )
    .set { ch_qc }
    ch_versions = ch_versions.mix(ch_qc.versions)

    SNV_CALLING (
        ch_mapped.bam_umi,
        ch_mapped.bam_dedup,
        beds,
        ch_qc.melt_qc,
        ch_qc.dedup_bam_is_metrics.groupTuple()
    )
    .set { ch_vcf }
    ch_versions = ch_versions.mix(ch_vcf.versions)

    SNV_ANNOTATE ( ch_vcf.agg_vcf, ch_vcf.concat_vcfs )
    .set { ch_vcf_anno }
    ch_versions = ch_versions.mix(ch_vcf_anno.versions)

    CNV_CALLING ( 
        ch_mapped.bam_umi, 
        ch_vcf_anno.germline_variants,
        ch_mapped.bam_dedup,
        gatk_ref
    )
    .set { ch_cnvcalled }
    ch_versions = ch_versions.mix(ch_cnvcalled.versions)

    BIOMARKERS ( 
        ch_cnvcalled.cnvkit_hrd,
        ch_mapped.bam_umi,
        ch_mapped.bam_dedup
    )
    .set { ch_bio }
    ch_versions = ch_versions.mix(ch_bio.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        CHECK_INPUT.out.meta
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
    logFile = file("${params.resultsdir}/cron/logs/" + base + ".complete")
    logFile.text = msg
    logFile.append(error)
}
