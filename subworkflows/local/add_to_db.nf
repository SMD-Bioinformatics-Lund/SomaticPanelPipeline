#!/usr/bin/env nextflow

include { COYOTE               } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf
        lowcov
        segments
        gens
        gatcov_plot

    main:

        COYOTE { vcf.join(lowcov).join(segments).join(gatcov_plot) }

    emit:
        coyotedone = COYOTE.out.coyote_import
        

}

