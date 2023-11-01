#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { CHECK_INPUT                     } from '../subworkflows/local/create_input_cnvrefs'
include { BED_INTERVALS                   } from '../subworkflows/local/bed_intervals'
include { CALL_COHORT                     } from '../subworkflows/local/call_cohort'
include { CNVKITREFS                      } from '../subworkflows/local/cnvkit_refs'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'

csv = file(params.csv)


workflow SPP_CREATE_CNVPON {

    ch_versions = Channel.empty()

    // Checks input, creates meta-channel and decides whether data should be downsampled //
    CHECK_INPUT ( Channel.fromPath(csv) )

    CNVKITREFS( 
        CHECK_INPUT.out.sample,
        params.name
    )
    ch_versions = ch_versions.mix(CNVKITREFS.out.versions)

    BED_INTERVALS(
        params.name,
        CHECK_INPUT.out.sample
    )
    .set { ch_bed }
    ch_versions = ch_versions.mix(ch_bed.versions)

    CALL_COHORT (
        ch_bed.intervals,
        ch_bed.intervals_scattered,
        ch_bed.counts
    )
    ch_versions = ch_versions.mix(CALL_COHORT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

