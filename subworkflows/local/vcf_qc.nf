#!/usr/bin/env nextflow

include { CONTAMINATION            } from '../../modules/local/qc/main'
include { PEDDY                    } from '../../modules/local/peddy/main'

workflow VCF_QC {
    take:        
        vep_vcf                    // channel: [ val(group), val(meta), file("*.vep.vcf") ]
        tumor_germline             // channel: [ val(group), file(vcf), file(tbi) ]
        meta_ch                    // channel: [mandatory] [ [sample_id, group, sex, phenotype, paternal_id, maternal_id, case_id] ]

    main:
        ch_versions = Channel.empty()

        vep_vcf_filtered = vep_vcf.filter { group, meta, vcf ->
            meta.sex != false
        }

        ped_ch = meta_ch.map { group, meta ->

            def sex_code = (
                meta.sex == 'M'   ? 1 :
                meta.sex == 'F' ? 2 : 0
            )

            def ped_line = "${group}\t${meta.id}\t0\t0\t${sex_code}\t2\n"

            def ped_file = file("${meta.id}.ped")
            ped_file.text = ped_line

            tuple(group, meta, ped_file)
        }

        ped_ch_tumor = ped_ch.filter { group, meta, ped ->
            meta.type == "T"
        }.join(tumor_germline).view()

        PEDDY { ped_ch_tumor }

        CONTAMINATION { vep_vcf }
        ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
        
    emit:
        qcdone                  =   CONTAMINATION.out.contamination_cdm                  // channel: [ tuple val(group), file("dist.txt"), file("sampleid.png") ]
        results                 =   CONTAMINATION.out.contamination_result_files         // channel: [ val(group), val(meta), file(cdm) ]

}