process UMI_CONSENSUS {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("consensus.fastq.gz"),    emit: consensus_fastq
        tuple val(group), val(meta), file("noumi.sam"),             emit: noumi_sam
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''                       
        def args2   = task.ext.args2    ?: ''
        def args3   = task.ext.args3    ?: ''
        """
        export skip_coord_end=true

        sentieon umi extract $args $r1 $r2 \\
        | sentieon bwa mem \\
            -t ${task.cpus} \\
            $args2 - \\
        | tee -a noumi.sam \\
        | sentieon umi consensus $args3 -o consensus.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch consensus.fastq.gz noumi.sam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}

process BWA_UMI {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(consensus_fastq)

    output:
        tuple val(group), val(meta), file("${out_umi}"), file("${out_umi}.bai"),    emit: bam_umi
        path "versions.yml",                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
        }

        """
        sentieon bwa mem \\
            -t ${task.cpus} \\
            $args \\
            ${consensus_fastq} \\
        | sentieon util sort -i - \\
            -o ${out_umi} \\
            $args2

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """

    stub:
        out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
        }

        """
        touch ${out_umi} ${out_umi}.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}

process BWA_SORT {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(sam)

    output:
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"),    emit: sorted_bam
        path "versions.yml",                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''

        out_bam = meta.id+"."+meta.type+".bwa.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
        }

        """
        sentieon util sort -i ${sam} -o ${out_bam} $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        out_bam = meta.id+"."+meta.type+".bwa.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
        }

        """
        touch ${out_bam} ${out_bam}.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process LOCUS_COLLECTOR {
    label 'process_very_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta),
              file(bam), file(bai),
              file("score.gz"), file("score.gz.tbi"), emit: score
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ""
        """
        sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        """
        touch score.gz score.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process MARKDUP {
    label 'process_very_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(score), file(score_tbi)

    output:
        tuple val(group), val(meta), file("*dedup_metrics.txt"),                 emit: metrics
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), emit: dedup_bam
        path "versions.yml",                                                     emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""
        def prefix  = task.ext.prefix   ?: ""

        out_bam = meta.id+"."+meta.type+".dedup.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
        }
        """
        sentieon driver -t ${task.cpus} -i $bam --algo Dedup $args --metrics ${prefix}dedup_metrics.txt $out_bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: ""

        out_bam = meta.id+"."+meta.type+".dedup.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
        }
        """
        touch ${out_bam} ${out_bam}.bai ${prefix}dedup_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process BQSR_UMI {
    label 'process_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file(bam), file(bai), file("*.bqsr.table"),   emit: bam_bqsr
        path "versions.yml",                                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference and common arguments for driver
        def args2   = task.ext.args2    ?: ""   // algo specific arguments
        def prefix  = task.ext.prefix   ?: "${meta.id}"

        """
        sentieon driver $args -t ${task.cpus} -i $bam --algo QualCal $args2 ${prefix}.bqsr.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}.bqsr.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process SENTIEON_QC {
    label 'process_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(dedup)

    output:
        tuple val(group), val(meta), file(bam), file(bai), file("*_is_metrics.txt"),               emit: dedup_bam_is_metrics
        tuple val(group), val(meta), file("*.txt"), file("cov_metrics.txt.sample_summary"),        emit: qc_files
        path "versions.yml",                                                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference 
        def args2   = task.ext.args2    ?: ""   // common arguments for driver
        def args3   = task.ext.args3    ?: ""
        def args4   = task.ext.args4    ?: ""
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        sentieon driver \\
            $args \\
            $args2 -t ${task.cpus} -i ${bam} \\
            --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \\
            --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt \\
            --algo InsertSizeMetricAlgo is_metrics.txt \\
            --algo CoverageMetrics $args3 cov_metrics.txt

        sentieon driver \\
            $args -t ${task.cpus} -i ${bam} \\
            --algo HsMetricAlgo $args4 hs_metrics.txt

        cp is_metrics.txt ${prefix}_is_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}_is_metrics.txt
        touch ${prefix}_${meta.type}.QC
        touch cov_metrics.txt.sample_summary

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process TNSCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        each file(bed)

    output:
        tuple val(group), val(meta), val("tnscope"), file("tnscope_${bed}.vcf.raw"),    emit: vcfparts_tnscope
        path "versions.yml",                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            sentieon driver -t ${task.cpus} $args \\
                -i ${bams[tumor_idx]} -q ${bqsr[tumor_idx]} \\
                -i ${bams[normal_idx]} -q ${bqsr[normal_idx]} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[tumor_idx]} --normal_sample ${meta.id[normal_idx]} \\
                $args2 \\
                --min_tumor_allele_frac ${params.tnscope_var_freq_cutoff_p} \\
                tnscope_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            sentieon driver -t ${task.cpus} $args \\
                -i ${bams} -q ${bqsr} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[0]} \\
                $args2 \\
                --min_tumor_allele_frac ${params.tnscope_var_freq_cutoff_up} \\
                tnscope_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """
        }

    stub:
        """
        touch tnscope_${bed}.vcf.raw

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process DNASCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: normal_germline
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        sentieon driver $args \\
            -i $bam \\
            --interval ${params.regions_bed} \\
            -q $bqsr \\
            --algo DNAscope \\
            ${meta.id}.dnascope.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.dnascope.vcf.gz
        touch ${meta.id}.dnascope.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}
