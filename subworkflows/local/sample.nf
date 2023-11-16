#!/usr/bin/env nextflow

include { SEQTK                } from '../../modules/local/seqtk/main'
include { FASTP                } from '../../modules/local/fastp/main'


workflow SAMPLE {
    take:
        fastq       // channel: [ val(meta), [ reads ] ]

    main:
        ch_versions = Channel.empty()

        // Only sub-sample if meta.sub == value
        SEQTK ( fastq.filter{ item -> item[1].sub != false } )
        ch_versions = ch_versions.mix(SEQTK.out.versions)


        //combine with any sample that does not get sub-sampled
        fastq_sample = SEQTK.out.fastq_sub.mix( fastq.filter{ item -> item[1].sub == false } )

        if (params.trimfq) {
            FASTP ( fastq_sample )
            fastq_done = FASTP.out.fastq_trimmed
            ch_versions = ch_versions.mix(FASTP.out.versions)
        }
        else {
            fastq_done = fastq_sample
        }

    emit:
        fastq_trim  =   fastq_done      // channel: [ val(meta), [ reads ] ] 
        versions    =   ch_versions     // channel: [ file(versions) ]

}