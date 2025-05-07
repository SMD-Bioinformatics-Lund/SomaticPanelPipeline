process BWA_UMI {
    label 'process_alot'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("${out_umi}"), file("${out_umi}.bai"),    emit: bam_umi
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"),    emit: bam_umi_markdup
        path "versions.yml",                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''                       
        def args2   = task.ext.args2    ?: ''
        def args3   = task.ext.args3    ?: ''
        def args4   = task.ext.args4    ?: ''
        def args5   = task.ext.args5    ?: ''

        out_bam = meta.id+"."+meta.type+".bwa.sort.bam"
        out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
            out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
        }

        """
        export skip_coord_end=true
        
        sentieon umi extract $args $r1 $r2 \\
        | sentieon bwa mem \\
            -t ${task.cpus} \\
            $args2 - \\
        | tee -a noumi.sam \\
        | sentieon umi consensus $args3 -o consensus.fastq.gz

        sentieon bwa mem \\
            -t ${task.cpus} \\
            $args4 \\
            consensus.fastq.gz \\
        | sentieon util sort -i - \\
            -o ${out_umi} \\
            $args5

        sentieon util sort -i noumi.sam -o ${out_bam} --sam2bam

        rm noumi.sam

        touch dedup_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """

    stub:
        out_bam = meta.id+"."+meta.type+".bwa.sort.bam"
        out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"

        if (meta.sub) {
            submbp  = params.sample_val / 1000000
            submbp  = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
            out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
        }

        """
        touch ${out_bam} ${out_bam}.bai
        touch ${out_umi} ${out_umi}.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}

process MARKDUP {
    label 'process_very_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*dedup_metrics.txt"),                                                emit: dedup_metrics
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"),                                emit: bam_qc
        path "versions.yml",                                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""                       
        def args2   = task.ext.args2    ?: ""
        def prefix  = task.ext.prefix   ?: ""

        out_bam = meta.id+"."+meta.type+".dedup.bam"

        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
        }
        """
        sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector $args
        sentieon driver -t ${task.cpus} -i $bam --algo Dedup $args2 --metrics ${prefix}dedup_metrics.txt $out_bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def args    = task.ext.args     ?: ""                       
        def args2   = task.ext.args2    ?: ""
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
        tuple val(group), val(meta), file(bam), file(bai), file("*.bqsr.table"),   emit: bam_varcall
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
        tuple val(group), val(meta), file(bam), file(bai), file("*_is_metrics.txt"),   emit: dedup_bam_is_metrics
        tuple val(group), val(meta), file("*.txt"),                                    emit: qc_files
        path "versions.yml",                                                           emit: versions

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process SENTIEON_QC_TO_CDM {
    label 'process_high'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(qc_files)

    output:
        tuple val(group), val(meta), file("*_${meta.type}.QC"),                        emit: qc_cdm

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        qc_sentieon.pl ${meta.id}_${meta.type} panel > ${prefix}_${meta.type}.QC
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}_${meta.type}.QC
        """
}

process TNSCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        each file(bed)

    output:
        tuple val(group), val(meta), file("tnscope_${bed}.vcf.raw"),        emit: vcfparts_tnscope
        path "versions.yml",                                                emit: versions

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
        touch tnscope_${bed}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process FILTER_TNSCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val("tnscope"), val(group), file("tnscope_${bed}.vcf"),   emit: vcfparts_tnscope_filtered
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            filter_tnscope_somatic.pl $vcf ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > tnscope_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            filter_tnscope_unpaired.pl $vcf > tnscope_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """ 
        }

    stub:
        """
        touch tnscope_${bed}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}