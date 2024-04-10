#!/usr/bin/env nextflow


include { ALLELE_CALL            } from '../../modules/local/idSnp/main'
include { SNP_CHECK              } from '../../modules/local/idSnp/main'


workflow ID_SNP {
    take:
        bam_dedup                                                   // channel: [ val(group), val(meta), file(bam), file(bai)] 


    main:
        ch_versions = Channel.empty()

        if( meta.id.size() >= 2 ) {
            ALLELE_CALL (bam_dedup)
            ch_versions = ch_versions.mix(ALLELE_CALL.out.versions)
            
            SNP_CHECK(ALLELE_CALL.out.sample_id_vcf.groupTuple())
            ch_versions = ch_versions.mix(SNP_CHECK.out.versions)
        }
        else {
            
        }
        

    emit:
        idsnp               =   SNP_CHECK.out.idsnp_checked         // channel: [ val(group), val(meta), file(snpid) ]
        versions            =   ch_versions                         // channel: [ file(versions) ]
}

