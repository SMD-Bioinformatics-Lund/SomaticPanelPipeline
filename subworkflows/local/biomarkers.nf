#!/usr/bin/env nextflow

include { ASCAT_252                   } from '../../modules/local/ASCAT/main'
include { ASCAT_30                    } from '../../modules/local/ASCAT/main'
include { CNVKIT2ASCAT                } from '../../modules/local/ASCAT/main'
include { MSISENSOR                   } from '../../modules/local/msisensor/main'
include { GENEFUSE                    } from '../../modules/local/genefuse/main'
include { CNVKIT2OVAHRDSCAR           } from '../../modules/local/hrdsw/main'
include { CNVKIT2SCARHRD              } from '../../modules/local/hrdsw/main'
include { ASCAT2SCARHRD               } from '../../modules/local/hrdsw/main'
include { ASCAT2OVAHRDSCAR            } from '../../modules/local/hrdsw/main'
include { SCARHRD                     } from '../../modules/local/hrdsw/main'
include { OVAHRDSCAR                  } from '../../modules/local/hrdsw/main'



workflow BIOMARKERS {
    take: 
        baflogr
        cnvkitsegments
        cnvkitsegment_purity

    main:
        CNVKIT2ASCAT ( baflogr )
        //ASCAT_252 ( CNVKIT2ASCAT.out.ascat_input )
        ASCAT_30 ( CNVKIT2ASCAT.out.ascat_input )
        CNVKIT2OVAHRDSCAR ( cnvkitsegments.mix(cnvkitsegment_purity) )
        CNVKIT2SCARHRD ( cnvkitsegments.mix(cnvkitsegment_purity) ) //.join(ASCAT_30.out.ploidy, by:[0,1])
        ASCAT2SCARHRD ( ASCAT_30.out.baflogr ) //.join(ASCAT_30.out.ploidy, by:[0,1])
        ASCAT2OVAHRDSCAR ( ASCAT_30.out.baflogr )
        SCARHRD ( ASCAT2SCARHRD.out.scarHRD_segments.mix(CNVKIT2SCARHRD.out.scarHRD_segments) )
        OVAHRDSCAR ( ASCAT2OVAHRDSCAR.out.ovaHRDscar_segments.mix(CNVKIT2OVAHRDSCAR.out.ovaHRDscar_segments) )

    emit:
        ascat = CNVKIT2ASCAT.out.ascat_input
        //hrdscore_ascat252 = ASCAT_252.out.ascat252_HRD_score
        //hrdscore_ascat30 = ASCAT_30.out.ascat30_HRD_score

}