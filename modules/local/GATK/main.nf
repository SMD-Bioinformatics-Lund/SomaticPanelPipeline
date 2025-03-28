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


process GATKCOV_COUNT {
    label 'process_high'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.standardizedCR.tsv"), file("*.denoisedCR.tsv"),    emit: gatk_count
        path "versions.yml",                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args      ?: ''
        def args2       = task.ext.args2     ?: ''
        def args3       = task.ext.args3     ?: ''
        def prefix      = task.ext.prefix    ?: "${meta.id}"
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
        gatk CollectReadCounts \\
            -I ${bam} \\
            $args \\
            -O ${bam}.hdf5

        gatk --java-options "-Xmx${avail_mem}M" DenoiseReadCounts \\
            -I ${bam}.hdf5 $args2 \\
            --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
            --denoised-copy-ratios ${prefix}.denoisedCR.tsv

        gatk PlotDenoisedCopyRatios \\
            --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
            --denoised-copy-ratios ${prefix}.denoisedCR.tsv \\
            $args3 \\
            --output . --output-prefix ${prefix}

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

        touch ${prefix}.standardizedCR.tsv 
        touch ${prefix}.denoisedCR.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """    
}


process GATKCOV_CALL {
    label 'process_high'
    tag "${meta.id[tumor_idx]}"

    input:
        tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

    output:
        tuple val(group), file("*.called.seg"),     emit: gatcov_called
        tuple val(group), file("*.modeled.png"),    emit: gatcov_plot
        path "versions.yml",                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args      ?: ''
        def args2  = task.ext.args2     ?: ''
        def args3  = task.ext.args3     ?: ''

        normalcounts = ""
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        if (meta.id.size() == 2) {
            args = args + " --normal-allelic-counts " + allele[normal_idx]
        }

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
        gatk --java-options "-Xmx${avail_mem}M" ModelSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${allele[tumor_idx]} \\
            $args \\
            --output . \\
            --output-prefix ${prefix}

        gatk CallCopyRatioSegments \\
            --input ${prefix}.cr.seg \\
            --output ${prefix}.called.seg \\
            $args2

        gatk PlotModeledSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${prefix}.hets.tsv \\
            --segments ${prefix}.modelFinal.seg \\
            $args3 \\
            --output . \\
            --output-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        prefix = task.ext.prefix ?: "${meta.id[tumor_idx]}"
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${prefix}.called.seg
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


process POSTPROCESS {
    label 'process_medium'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(i), file(tar), file(ploidy), val(shard_no), val(shard)

    output:
        tuple val(group), val(meta), file("*_gatk_genotyped-intervals.vcf.gz"),    emit: gatk_germline_intervals
        tuple val(group), val(meta), file("*_gatk_genotyped-segments.vcf.gz"),     emit: gatk_germline_segmentsvcf
        tuple val(group), val(meta), file("*_gatk_denoised.vcf.gz"),               emit: gatk_germline_denoised_log2
        path "versions.yml",                                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        args        = task.ext.args     ?: ''                   
        prefix      = task.ext.prefix   ?: "${meta.id}"
        avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        modelshards = shard.join(' --model-shard-path ') // join each reference shard
        caseshards = []
        for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
            tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls' 
            caseshards = caseshards + tmp
        }
        caseshards = caseshards.join( ' --calls-shard-path ')

    shell:
        '''
        export THEANO_FLAGS="base_compiledir=."
        for model in !{tar}; do
        tar -xvf $model
        done
        tar -xvf !{ploidy}
        set +u
        source activate gatk
        export MKL_NUM_THREADS=!{task.cpus}
        export OMP_NUM_THREADS=!{task.cpus}
        gatk --java-options "-Xmx!{avail_mem}M" PostprocessGermlineCNVCalls \
            --contig-ploidy-calls ploidy/!{group}-calls/ \
            --output-genotyped-intervals !{prefix}_gatk_genotyped-intervals.vcf.gz \
            --output-genotyped-segments !{prefix}_gatk_genotyped-segments.vcf.gz \
            --output-denoised-copy-ratios !{prefix}_gatk_denoised.vcf.gz \
            !{args} \
            --calls-shard-path !{caseshards} \
            --model-shard-path !{modelshards}

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            gatk4: $(echo $(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*$//')
            tar: $(echo $(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        '''

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
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process FILTER_MERGE_GATK {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), file("*.gatk.filtered.merged.vcf"),   emit: gatk_normal_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        filter_gatk.pl $gatk > ${prefix}.gatk.filtered.vcf
        mergeGATK.pl ${prefix}.gatk.filtered.vcf > ${prefix}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}


process GATK2VCF {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), file(seg), val(meta)

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
        mergeGATK_tumor.pl $gatk > ${prefix}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}