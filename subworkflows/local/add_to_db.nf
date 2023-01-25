#!/usr/bin/env nextflow

include { COYOTE               } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf
        lowcov
        segments
        gens

    main:

        COYOTE { vcf.join(lowcov).join(segments) }

    emit:
        coyotedone = COYOTE.out.coyote_import
        

}

