#!/usr/bin/env nextflow

include { COYOTE               } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf             // val(group), val(meta), file(vcf)
        lowcov          // val(group), val(meta.type), file(lowcov)
        segments        // val(group), file(segments)
        gens            // val(group), val(meta), file(gens)
        gatcov_plot     // val(group), file(plot)
        fusions         // val(group), file(vcf)
        biomarkers      // val(group), file(json)

    main:
        lc = lowcov.map{ val-> tuple(val[0], val[2] ) }
        optional = lc.mix(segments,gatcov_plot,biomarkers,fusions).groupTuple()
        COYOTE { vcf.join(optional) }

    emit:
        coyotedone = COYOTE.out.coyote_import
        

}

