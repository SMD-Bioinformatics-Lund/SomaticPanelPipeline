#!/usr/bin/env nextflow

include { COYOTE               } from '../../modules/local/coyote/main'

workflow ADD_TO_DB {
    take: 
        vcf             // channel: [mandatory] [ val(group), val(meta), file(vcf) ]
        lowcov          // channel: [mandatory] [ val(group), val(meta.type), file(lowcov) ]
        segments        // channel: [optional] [ val(group), file(segments) ]
        gens            // channel: [optional] [ val(group), val(meta), file(gens) ]
        gatcov_plot     // channel: [optional] [ val(group), file(plot) ]
        fusions         // channel: [optional] [ val(group), file(vcf) ]
        biomarkers      // channel: [optional] [ val(group), file(json) ]

    main:
        lc = lowcov.map{ val-> tuple(val[0], val[2] ) }
        optional = lc.mix(segments,gatcov_plot,biomarkers,fusions).groupTuple()
        COYOTE { vcf.join(optional) }

    emit:
        coyotedone = COYOTE.out.coyote_import   // channel: [ val(group), file(coyote) ]
        
}

