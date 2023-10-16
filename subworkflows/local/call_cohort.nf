#!/usr/bin/env nextflow

include { COHORT_PLOIDY           } from '../../modules/local/GATK_references/main'
include { COHORT_CALL             } from '../../modules/local/GATK_references/main'
include { COHORT_CALL_PANEL       } from '../../modules/local/GATK_references/main'
include { GATK_SOM_PON            } from '../../modules/local/GATK_references/main'


workflow CALL_COHORT {
    take:
        intervals           // val(prefix), file(intervals) name of reference
        intervals_scattered // val(prefix), file(scatters.tar)
        counts              // val(prefix), file(id), file(tsv)

    main:
        ch_versions = Channel.empty()

        COHORT_PLOIDY(intervals,counts.groupTuple())
        scatters    = intervals_scattered.map{ val-> val[1] }.flatten()
        ch_versions = ch_versions.mix(COHORT_PLOIDY.out.versions)

        if (params.panel) {
            COHORT_CALL_PANEL(counts.groupTuple().join(COHORT_PLOIDY.out.ploidy).join(intervals))
            ch_versions = ch_versions.mix(COHORT_CALL_PANEL.out.versions)
        }
        else {
            COHORT_CALL(counts.groupTuple().join(COHORT_PLOIDY.out.ploidy).combine(scatters))
            ch_versions = ch_versions.mix(COHORT_CALL.out.versions)
        }

        GATK_SOM_PON( counts.groupTuple() )
        ch_versions = ch_versions.mix(GATK_SOM_PON.out.versions)

    emit:
        ploidy      =   COHORT_PLOIDY.out.ploidy
        versions    =   ch_versions
}