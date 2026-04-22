#!/usr/bin/env nextflow

include { CONTAMINATION            } from '../../modules/local/qc/main'
include { PEDDY                    } from '../../modules/local/peddy/main'

workflow VCF_QC {
    take:        
        vep_vcf                    // channel: [ val(group), val(meta), file("*.vep.vcf") ]

    main:
        ch_versions = Channel.empty()

        vep_vcf_filtered = vep_vcf.filter { group, meta, vcf ->
            meta.sex != false
        }

        ped_ch = vep_vcf.map { group, meta, vcf ->

            def sex_value = meta.containsKey('sex') ? meta.sex : null

            def sex_code = (
                sex_value == 'M'   ? 1 :
                sex_value == 'F' ? 2 : 0
            )

            def ped_line = "${group} ${meta.id} 0 0 ${sex_code} 2\n"

            def ped_file = file("${meta.id}.ped")
            ped_file.text = ped_line

            tuple(group, meta, ped_file)
        }


        PEDDY { vep_vcf_filtered.join(ped_ch, by:[0,1]) }

        CONTAMINATION { vep_vcf }
        ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
        
    emit:
        qcdone                  =   CONTAMINATION.out.contamination_cdm                  // channel: [ tuple val(group), file("dist.txt"), file("sampleid.png") ]
        results                 =   CONTAMINATION.out.contamination_result_files         // channel: [ val(group), val(meta), file(cdm) ]

}