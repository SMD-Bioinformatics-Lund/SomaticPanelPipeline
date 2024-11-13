process PREPROCESSINTERVALS {
    label 'process_low'

    input:
        val(reference)

    output:
        tuple val(reference), file("${reference}.preprocessed.blacklisted.interval_list"),    emit: preprocessed
        path "versions.yml",                                                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def panel   = params.panel ? "-L ${params.interval_list}" : ''

        """
        gatk PreprocessIntervals $args -O ${reference}.preprocessed.blacklisted.interval_list $panel

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch ${reference}.preprocessed.blacklisted.interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process COUNT_READS {
    label 'process_low'

    input:
        tuple val(id), file(cram), file(crai), file(bai)
        tuple val(reference), path(bed)

    output:
        tuple val(reference), val(id), file("*.tsv"),  emit: count_tsv
        path "versions.yml",                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${id}"
        """
        gatk CollectReadCounts $args -O ${prefix}.tsv -I $cram -L $bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${id}"
        """
        touch ${prefix}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process ANNOTATE_GC {
    label 'process_low'

    input:
        tuple val(reference), path(interval_list)

    output:
        tuple val(reference), file("*.annotated.tsv"), emit: annotated_intervals
        path "versions.yml",                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"
        """
        gatk AnnotateIntervals -L $interval_list $args -O ${prefix}.annotated.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"
        """
        touch ${prefix}.annotated.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process CORRECT_GC {
    label 'process_min_cpus'
    label 'process_very_high_memory'

    input:
        tuple val(reference), file(interval_list), file(annotated)
        tuple val(reference), val(id), file(tsvs)

    output:
        tuple val(reference), file("*.preprocessed.blacklisted.gcfiltered.interval_list"), emit: corrected_intervals
        path "versions.yml",                                                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        gatk FilterIntervals -L $interval_list \\
            --annotated-intervals $annotated \\
            $args \\
            -O ${prefix}.preprocessed.blacklisted.gcfiltered.interval_list \\
            $tsv_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        touch ${prefix}.preprocessed.blacklisted.gcfiltered.interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process REMOVE_ALT_CONTIGS {
    label 'process_low'

    input:
        tuple val(reference), path(interval_list)

    output:
        tuple val(reference), path("*.preprocessed.blacklisted.gcfiltered.noalt.interval_list"), emit: noaltcontigs

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${reference}"
        """
        grep -vP "^[A-z,0-9]+_" $interval_list > ${prefix}.preprocessed.blacklisted.gcfiltered.noalt.interval_list
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"
        """
        touch ${prefix}.preprocessed.blacklisted.gcfiltered.noalt.interval_list
        """
}

process SCATTER_INTERVALS {
    label 'process_low'

    input:
        tuple val(reference), path(interval_list)
        val(size)

    output:
        tuple val(reference), path("scatter/*"),    emit: scatter
        path "versions.yml",                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        count = 0
        """
        mkdir scatter
        gatk IntervalListTools $args --INPUT $interval_list --SCATTER_CONTENT $size --OUTPUT scatter
        ls scatter/*/scattered.interval_list | wc -l

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir scatter
        mkdir scatter/temp_001
        mkdir scatter/temp_002
        touch scatter/temp_001/scattered.interval_list
        touch scatter/temp_002/scattered.interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process COHORT_PLOIDY {
    label 'process_medium'

    input:
        tuple val(reference), file(interval_list)
        tuple val(reference), val(id), file(tsvs)

    output:
        tuple val(reference), path("*_ploidy-model"), path("*_ploidy-calls"),  emit: ploidy
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        gatk --java-options "-Djava.io.tmpdir=/local/scratch/" \\
            DetermineGermlineContigPloidy \\
            -L $interval_list \\
            $args \\
            --output . \\
            --output-prefix ${prefix}_ploidy \\
            $tsv_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        mkdir ${prefix}_ploidy-model ${prefix}_ploidy-calls

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process COHORT_CALL {
    label 'process_medium_cpus'
    label 'process_very_high_memory'

    input:
        tuple val(reference), val(id), file(tsvs), path(ploidy_model), path(ploidy_calls), file(scatter)

    output:
        tuple val(reference), path("cohort_calls/cohort_${name}"),  emit: calls
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        name = scatter.toString()
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        gatk --java-options "-Djava.io.tmpdir=/local/" \\
            GermlineCNVCaller --run-mode COHORT \\
            -L ${scatter}/scattered.interval_list \\
            --contig-ploidy-calls $ploidy_calls \\
            $args \\
            --output cohort_calls \\
            --output-prefix cohort_${name} \\
            $tsv_list
        touch test

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        name = scatter.toString()
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        echo $tsv_list
        mkdir cohort_calls
        mkdir cohort_calls/cohort_${name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process COHORT_CALL_PANEL {
    label 'process_medium_cpus'
    label 'process_very_high_memory'

    input:
        tuple val(reference), val(id), file(tsvs), path(ploidy_model), path(ploidy_calls), file(intervals)

    output:
        tuple val(reference), path("cohort_calls/${reference}*"),  emit: calls
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        gatk --java-options "-Djava.io.tmpdir=/local/scratch/" \\
            GermlineCNVCaller --run-mode COHORT \\
            -L $intervals \\
            --contig-ploidy-calls $ploidy_calls \\
            $args \\
            --output cohort_calls \\
            --output-prefix ${reference} \\
            $tsv_list
        touch test

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        echo $tsv_list
        mkdir cohort_calls
        mkdir cohort_calls/${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

}

process GATK_SOM_PON {
    label 'process_medium_cpus'
    label 'process_very_high_memory'

    input:
        tuple val(reference), val(id), file(tsvs)

    output:
        tuple val(reference), path("*.somatic_gatk_pon.hdf5"), emit: somatic_pon
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        gatk --java-options "-Djava.io.tmpdir=/local/scratch/" \\
            CreateReadCountPanelOfNormals \\
            $args \\
            -O ${prefix}.somatic_gatk_pon.hdf5 \\
            $tsv_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${reference}"

        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        echo $tsv_list
        touch ${prefix}.somatic_gatk_pon.hdf5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}