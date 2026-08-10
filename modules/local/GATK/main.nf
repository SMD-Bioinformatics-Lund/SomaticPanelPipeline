// SOMATIC CALLING //
process GATKCOV_BAF {
    label 'process_medium_cpus'
    label 'process_high_memory'
    label 'process_more_time'

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.allelicCounts.tsv"),  emit: gatk_baf
        path "versions.yml",                                       emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.id}"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK CollectAllelicCounts] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" CollectAllelicCounts \\
            -I $bam \\
            $args \\
            -O ${prefix}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATK_COLLECT_READ_COUNTS {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.hdf5"), emit: read_counts
        path "versions.yml",                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk CollectReadCounts \\
            -I ${bam} \\
            $args \\
            -O ${prefix}.hdf5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.hdf5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK_DENOISE_READ_COUNTS {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(hdf5)

    output:
        tuple val(group), val(meta), file("*.standardizedCR.tsv"), file("*.denoisedCR.tsv"), emit: denoised_counts
        path "versions.yml",                                                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args      = task.ext.args   ?: ''
        def prefix    = task.ext.prefix ?: "${meta.id}"
        def avail_mem = 12288
        if (task.memory) {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk --java-options "-Xmx${avail_mem}M" DenoiseReadCounts \\
            -I ${hdf5} $args \\
            --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
            --denoised-copy-ratios ${prefix}.denoisedCR.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.standardizedCR.tsv
        touch ${prefix}.denoisedCR.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK_PLOT_DENOISED_COPY_RATIOS {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(stdCR), file(denoised)

    output:
        tuple val(group), val(meta), file("*.denoised.png"), optional: true, emit: denoised_plot
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk PlotDenoisedCopyRatios \\
            --standardized-copy-ratios ${stdCR} \\
            --denoised-copy-ratios ${denoised} \\
            $args \\
            --output . --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.denoised.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATK_MODEL_SEGMENTS {
    label 'process_high'
    tag "$group"

    input:
        tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

    output:
        tuple val(group), val(meta), file("*.cr.seg"), file("*.hets.tsv"), file("*.modelFinal.seg"), emit: model_segments
        path "versions.yml",                                                                                 emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        if (meta.id.size() == 2) {
            args = args + " --normal-allelic-counts " + allele[normal_idx]
        }
        def avail_mem = 12288
        if (task.memory) {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk --java-options "-Xmx${avail_mem}M" ModelSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${allele[tumor_idx]} \\
            $args \\
            --output . \\
            --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        touch ${prefix}.cr.seg
        touch ${prefix}.hets.tsv
        touch ${prefix}.modelFinal.seg

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK_CALL_COPY_RATIO_SEGMENTS {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cr_seg), file(hets), file(model_final)

    output:
        tuple val(group), val(meta), file("*.called.seg"), emit: called_segments
        path "versions.yml",                              emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk CallCopyRatioSegments \\
            --input ${cr_seg} \\
            --output ${prefix}.called.seg \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.called.seg

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK_PLOT_MODELED_SEGMENTS {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(denoised), file(cr_seg), file(hets), file(model_final)

    output:
        tuple val(group), val(meta), file("*.modeled.png"), emit: modeled_plot
        path "versions.yml",                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        gatk PlotModeledSegments \\
            --denoised-copy-ratios ${denoised} \\
            --allelic-counts ${hets} \\
            --segments ${model_final} \\
            $args \\
            --output . \\
            --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.modeled.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


// GERMLINE CALLING //
process GATK_COUNT_GERMLINE {
    label 'process_high'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.tsv"),   emit: count_germline
        path "versions.yml",                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args                  ?: ''
        def prefix      = task.ext.prefix                ?: "${meta.id}"
        def suffix      = args.contains("--format TSV")  ? '.tsv' : '.tsv'
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" CollectReadCounts \\
            -I $bam \\
            $args \\
            -O ${prefix}${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args   = task.ext.args                  ?: ''
        def prefix = task.ext.prefix                ?: "${meta.id}"
        def suffix = args.contains("--format TSV")  ? '.tsv' : '.tsv'
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATK_CALL_PLOIDY {
    label 'process_high'
    label "stage"
    label "scratch"
    tag "$group"

    input:
        tuple val(group), val(meta), file(tsv)

    output:
        tuple val(group), val(meta), file("ploidy.tar"),    emit: gatk_ploidy
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args      ?: ''
        def prefix      = task.ext.prefix    ?: "${group}"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx${avail_mem}M" DetermineGermlineContigPloidy \\
            $args \\
            -I $tsv \\
            -O ploidy/ \\
            --output-prefix $prefix
        tar -cvf ploidy.tar ploidy/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix    ?: "${group}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ploidy.tar

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process GATK_CALL_GERMLINE_CNV {
    label 'process_medium'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(tsv), file(ploidy), val(i), val(refpart)

    output:
        tuple val(group), val(meta), val(i), file("*_${i}.tar"),    emit: gatk_call_germline
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args      ?: ''
        def prefix      = task.ext.prefix    ?: "${group}"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        tar -xvf ploidy.tar
        mkdir ${prefix}_${i}
        gatk --java-options "-Xmx${avail_mem}M" GermlineCNVCaller \\
            $args \\
            -I $tsv \\
            --contig-ploidy-calls ploidy/${prefix}-calls/ \\
            --model ${refpart} \\
            --output ${prefix}_${i}/ \\
            --output-prefix ${prefix}_${i}
        tar -cvf ${prefix}_${i}.tar ${prefix}_${i}/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix    ?: "${group}"
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        touch ${prefix}_${i}.tar

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process GATK_EXTRACT_GERMLINE_CNV_SHARDS {
    label 'process_low'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(i), file(tar), file(ploidy), val(shard_no), val(shard)

    output:
        tuple val(group), val(meta), val(i), path("calls"), path("ploidy"), val(shard_no), val(shard), emit: extracted_germline_cnv
        path "versions.yml",                                                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        extract_gatk_cnv_shards.sh \\
            --ploidy ${ploidy} \\
            --calls-dir calls \\
            ${tar}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p calls ploidy

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process POSTPROCESS {
    label 'process_medium'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(i), path(calls), path(ploidy), val(shard_no), val(shard)

    output:
        tuple val(group), val(meta), file("*_gatk_genotyped-intervals.vcf.gz"),    emit: gatk_germline_intervals
        tuple val(group), val(meta), file("*_gatk_genotyped-segments.vcf.gz"),     emit: gatk_germline_segmentsvcf
        tuple val(group), val(meta), file("*_gatk_denoised.vcf.gz"),               emit: gatk_germline_denoised_log2
        path "versions.yml",                                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args     ?: ''
        def prefix      = task.ext.prefix   ?: "${meta.id}"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        def modelshards = shard.join(' --model-shard-path ')
        def caseshards  = i.collect { shard_id -> "calls/${group}_${shard_id}/${group}_${shard_id}-calls" }.join(' --calls-shard-path ')

        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        gatk --java-options "-Xmx${avail_mem}M" PostprocessGermlineCNVCalls \\
            --contig-ploidy-calls ploidy/${group}-calls/ \\
            --output-genotyped-intervals ${prefix}_gatk_genotyped-intervals.vcf.gz \\
            --output-genotyped-segments ${prefix}_gatk_genotyped-segments.vcf.gz \\
            --output-denoised-copy-ratios ${prefix}_gatk_denoised.vcf.gz \\
            ${args} \\
            --calls-shard-path ${caseshards} \\
            --model-shard-path ${modelshards}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix    ?: "${meta.id}"

        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}_gatk_genotyped-intervals.vcf.gz 
        touch ${prefix}_gatk_genotyped-segments.vcf.gz 
        touch ${prefix}_gatk_denoised.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process FILTER_GATK {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), val(meta), file("*.gatk.filtered.vcf"), emit: gatk_filtered_vcf
        path "versions.yml",                                     emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        filter_gatk.py --vcf $gatk --out ${prefix}.gatk.filtered.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.gatk.filtered.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process MERGE_GATK {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), val(meta), file("*.gatk.filtered.merged.vcf"), emit: gatk_normal_vcf
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mergeGATK.py --vcf $gatk --out ${prefix}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process GATK2VCF {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(seg)

    output:
        tuple val(group), val(meta), file("*_gatk_tumor.vcf"),  emit: tumor_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def prefix = task.ext.prefix    ?: "${meta.id}"

        """
        gatk_to_vcf.py -s $seg -o ${prefix}_gatk_tumor.vcf $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_gatk_tumor.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process MERGE_GATK_TUMOR {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), val(meta), file("*_gatk_tumor_merged.vcf"),   emit: tumor_vcf_merged
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mergeGATK_tumor.py --vcf $gatk --out ${prefix}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}
