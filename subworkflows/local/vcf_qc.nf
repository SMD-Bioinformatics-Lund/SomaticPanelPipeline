#!/usr/bin/env nextflow

include { CONTAMINATION            } from '../../modules/local/qc/main'

workflow QC {
    take:        
        vep_vcf                    // channel: [ val(group), val(meta), file("*.vep.vcf") ]

    main:
        ch_versions = Channel.empty()

        CONTAMINATION { vep_vcf }
        ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
        
    emit:
        qcdone                  =   CONTAMINATION.out.contamination_cdm                  // channel: [ tuple val(group), file("dist.txt"), file("sampleid.png") ]
        results                 =   CONTAMINATION.out.contamination_result_files         // channel: [ val(group), val(meta), file(cdm) ]

}