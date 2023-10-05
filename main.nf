#!/usr/bin/env nextflow

/*
    A common somatic panel pipeline (SPP) for targeted gene panel analysis of myeloid, lymphoid and solid cancer data. 
*/

nextflow.enable.dsl = 2

include { SPP_CREATE_SNVPON        } from './workflows/pon_solid.nf'
include { SPP_COMMON               } from './workflows/commom.nf'
include { SPP_CREATE_CNVPON        } from './workflows/create_cnv_pons.nf'

println(params.genome_file)
csv = file(params.csv)
println(csv)

/*
    Different workflows for the analysis; current entry point is SPP i.e -entry SPP
*/
workflow SPP {
    SPP_COMMON ()
}

workflow SNVPON {
    SPP_CREATE_SNVPON()
}

workflow CNVPON {
    SPP_CREATE_CNVPON()
}


