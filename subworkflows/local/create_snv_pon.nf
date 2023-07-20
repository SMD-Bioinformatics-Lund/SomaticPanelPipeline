#!/usr/bin/env nextflow

include { CREATE_SNVPON               } from '../../modules/local/filters/main'

workflow CREATE_SNV_PON {
    take: 
        concat_vcfs

    main:
        CREATE_SNVPON( concat_vcfs.groupTuple(by:[1]))
    

}