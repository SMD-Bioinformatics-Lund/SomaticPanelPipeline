#!/usr/bin/env nextflow


nextflow.enable.dsl = 2


include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { SNV_CALLING                   } from '../subworkflows/local/snv_calling'
include { CNV_CALLING                   } from '../subworkflows/local/cnv_calling'
include { BIOMARKERS                    } from '../subworkflows/local/biomarkers'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

println(params.genome_file)
genome_file = file(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

csv = file(params.csv)
println(csv)

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, file(row.read1), file(row.read2)) }
    .set { fastq }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, (row.containsKey("ffpe") ? row.ffpe : false)) }
    .set { meta }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, (row.containsKey("purity") ? row.purity : false)) }
    .set { meta_purity }


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

    ALIGN_SENTIEON ( fastq )
    .set { ch_mapped }
    ch_versions = ch_versions.mix(ch_mapped.versions)

    SNV_CALLING ( ch_mapped.bam_umi.groupTuple(), beds, meta )
    .set { ch_vcf }
    ch_versions = ch_versions.mix(ch_vcf.versions)

    CNV_CALLING ( 
        ch_mapped.bam_umi, 
        ch_vcf.germline_variants,
        meta_purity,
        ch_vcf.concat_vcfs
    )
    .set { ch_cnvcalled }
    ch_versions = ch_versions.mix(ch_cnvcalled.versions)

    BIOMARKERS ( 
        ch_cnvcalled.baflogr,
        ch_cnvcalled.cnvkitsegment,
        ch_cnvcalled.cnvkitsegment_purity
    )
    ch_versions = ch_versions.mix(ch_bio.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
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
