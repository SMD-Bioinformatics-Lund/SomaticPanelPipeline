#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { SPP_SNVPON        } from './workflows/pon_solid.nf'
include { SOLID_GMS         } from './workflows/test_checkinput.nf'
include { CREATE_REF        } from './workflows/create_cnv_pons.nf'

println(params.genome_file)
csv = file(params.csv)
println(csv)

workflow SPP {
    SOLID_GMS()
}

workflow SNVPON {
    SPP_SNVPON()
}

workflow CNVPON {
    CREATE_REF()
}
