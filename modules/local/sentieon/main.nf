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
        params.umi

    script:
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
        
        sentieon umi extract -d 3M2S+T,3M2S+T $r1 $r2 \\
        |sentieon bwa mem \\
            -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -p -C $params.genome_file - \\
        |tee -a noumi.sam \\
        |sentieon umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

        sentieon bwa mem \\
            -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:illumina" \\
            -t ${task.cpus} \\
            -p -C $params.genome_file consensus.fastq.gz \\
        |sentieon util sort -i - \\
            -o ${out_umi} \\
            --sam2bam --umi_post_process

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
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
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
    tag "${meta.id}"
    label 'process_high'
    label 'scratch'
    label 'stage'
    
    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"),                            emit: bam_bqsr
        tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), file("dedup_metrics.txt"), emit: bam_qc
        path "versions.yml",                                                                                emit: versions

    script:
        out_bam = meta.id+"."+meta.type+".dedup.bam"
        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
        }
        """
        sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info score.gz
        sentieon driver -t ${task.cpus} -i $bam --algo Dedup --score_info score.gz --metrics dedup_metrics.txt $out_bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        out_bam = meta.id+"."+meta.type+".dedup.bam"
        if (meta.sub) {
            submbp = params.sample_val / 1000000
            submbp = submbp + "M"
            out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
        }
        """
        touch ${out_bam} ${out_bam}.bai dedup_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process BQSR_UMI {
    label 'process_medium'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"


    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file(bam), file(bai), file("${meta.id}.bqsr.table"),   emit: bam_varcall
        path "versions.yml",                                                                emit: versions

    when:
        params.umi

    script:
        """
        sentieon driver -t ${task.cpus} -r $params.genome_file -i $bam --algo QualCal ${meta.id}.bqsr.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.bqsr.table

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
        tuple val(group), val(meta), file(bam), file(bai), file("${meta.id}_is_metrics.txt"),   emit: dedup_bam
        tuple val(group), val(meta), file("${meta.id}_${meta.type}.QC"),                        emit: qc_cdm
        tuple val(group), val(meta), file("${meta.id}_${meta.type}.QC"),                        emit: qc_melt
        path "*.txt",                                                                           emit: txt
        path "versions.yml",                                                                    emit: versions

    script:
        """
        sentieon driver \\
            --interval $params.regions_bed_qc -r $params.genome_file -t ${task.cpus} -i ${bam} \\
            --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \\
            --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt \\
            --algo InsertSizeMetricAlgo is_metrics.txt \\
            --algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
        sentieon driver \\
            -r $params.genome_file -t ${task.cpus} -i ${bam} \\
            --algo HsMetricAlgo --targets_list $params.interval_list_qc --baits_list $params.interval_list_qc hs_metrics.txt

        cp is_metrics.txt ${meta.id}_is_metrics.txt

        qc_sentieon.pl ${meta.id}_${meta.type} panel > ${meta.id}_${meta.type}.QC

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_is_metrics.txt
        touch ${meta.id}_${meta.type}.QC

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}

process TNSCOPE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        each file(bed)

    output:
        tuple val("tnscope"), val(group), file("tnscope_${bed}.vcf"),   emit: vcfparts_tnscope
        path "versions.yml",                                            emit: versions

    when:
        params.tnscope

    script:
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            sentieon driver -t ${task.cpus} \\
                -r $params.genome_file \\
                -i ${bams[tumor_idx]} -q ${bqsr[tumor_idx]} \\
                -i ${bams[normal_idx]} -q ${bqsr[normal_idx]} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[tumor_idx]} --normal_sample ${meta.id[normal_idx]} \\
                --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
                --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
                tnscope_${bed}.vcf.raw

            filter_tnscope_somatic.pl tnscope_${bed}.vcf.raw ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > tnscope_${bed}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            sentieon driver -t ${task.cpus} -r $params.genome_file \\
                -i ${bams} -q ${bqsr} \\
                --interval $bed --algo TNscope \\
                --tumor_sample ${meta.id[0]} \\
                --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
                --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
                tnscope_${bed}.vcf.raw

            filter_tnscope_unpaired.pl tnscope_${bed}.vcf.raw > tnscope_${bed}.vcf

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