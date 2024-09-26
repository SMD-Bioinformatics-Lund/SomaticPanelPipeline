#!/usr/bin/env nextflow

include { ALLELE_CALL            } from '../../modules/local/idSnp/main'
include { SNP_CHECK              } from '../../modules/local/idSnp/main'
include { PAIR_CDM               } from '../../modules/local/idSnp/main'

workflow ID_SNP {
    take:
        bam_dedup    // channel: [ val(group), val(meta), file(bam), file(bai),val(pairedness)] 
        meta        // channel: [ id, group, phenotype, paternal_id, maternal_id, case_id ]                                        
        
    main:
        ch_versions = Channel.empty()

        ALLELE_CALL (bam_dedup)
        ch_versions = ch_versions.mix(ALLELE_CALL.out.versions)

        SNP_CHECK(ALLELE_CALL.out.sample_id_vcf.groupTuple())
        ch_versions = ch_versions.mix(SNP_CHECK.out.versions)

        PAIR_CDM (SNP_CHECK.out.idsnp_checked)

    emit:
        genotype            =   ALLELE_CALL.out.sample_id_vcf       // channel: [ val(group), val(meta), file("*.genotypes") ]
        idsnp               =   SNP_CHECK.out.idsnp_checked         // channel: [ val(group), val(meta), file(snpid) ]
        versions            =   ch_versions                         // channel: [ file(versions) ]
}

