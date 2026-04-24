#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
    take:
        csv
        paired

    main:
        CSV_CHECK(csv)

        checkedCSV = CSV_CHECK.out.csv
            .splitCsv(header:true, sep:',')
            .set { csvmap }

        // build meta once
        meta_ch = csvmap.map { row ->
            def meta = build_meta(row, paired)
            tuple(row.group, meta, row)
        }

        // separate pure meta channel
        meta = meta_ch.map { group, meta, row ->
            tuple(group, meta)
        }

        // build reads from same meta
        reads = meta_ch.map { group, meta_data, row ->
            tuple(group, meta, file(row.read1), file(row.read2))
        }

        fastq = reads.filter { group, meta_data, r1, r2 ->
            r1.name.endsWith("fastq.gz") || r1.name.endsWith("fq.gz")
        }

        bam = reads.filter { group, meta_data, r1, r2 ->
            r1.name.endsWith("bam") &&
            (r2.name.endsWith("bai") || r2.name.endsWith("bam.bai"))
        }

        vcf = reads.filter { group, meta_data, r1, r2 ->
            def f1 = r1.name.toLowerCase()
            def f2 = r2.name.toLowerCase()
            f1.endsWith("vcf") &&
            (f2.endsWith("tbi") || f2.endsWith("csi") || f2.endsWith("vcf.gz.tbi"))
        }

    emit:
        fastq
        bam
        vcf
        meta
}


def build_meta(LinkedHashMap row, paired) {
    def meta = [:]
    meta.id                 = row.id
    meta.group              = row.group
    meta.diagnosis          = row.diagnosis
    meta.type               = row.type
    meta.clarity_sample_id  = row.clarity_sample_id
    meta.ffpe               = row.containsKey("ffpe") && row.ffpe ? true : false
    meta.purity             = row.containsKey("purity") ? row.purity : false
    meta.sequencing_run     = row.sequencing_run
    meta.reads              = row.containsKey("n_reads") ? row.n_reads : false
    meta.sex                = row.containsKey("sex") ? row.sex : false
    meta.clarity_pool_id    = row.clarity_pool_id
    meta.paired             = paired

    // compute sub ONCE
    def sub = false
    if (meta.reads && params.sample) {  
        if (meta.reads.toInteger() > params.sample_val) {
            sub = (params.sample_val / meta.reads.toInteger()).round(2)
            if (sub == 1.00) sub = 0.99
        }
    }
    meta.sub = sub

    return meta
}
