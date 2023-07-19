process PREPROCESSINTERVALS {
    label 'gatk_small'

    input:
        val(prefix)

    output:
        tuple val(prefix), file("${prefix}.preprocessed.blacklisted.interval_list"), emit: preprocessed

    script:
        panel = ""
        if (params.panel) {
            panel = "-L ${params.interval_list}"
        }
        """
        gatk PreprocessIntervals -R ${params.genome_file} --padding ${params.padding} -imr OVERLAPPING_ONLY -O ${prefix}.preprocessed.blacklisted.interval_list -XL ${params.blacklist} $panel
        """

    stub:
        """
        touch ${prefix}.preprocessed.blacklisted.interval_list
        """
}

process COUNT_READS {
    
    input:
        tuple val(id), file(cram), file(crai), file(bai)
        tuple val(prefix), path(bed)

    output:
        tuple val(prefix), val(id), file("${id}.tsv"), emit: count_tsv

    script:
        """
        gatk CollectReadCounts -R ${params.genome_file} -imr OVERLAPPING_ONLY --format TSV -O ${id}.tsv -I $cram -L $bed
        """

    stub:
        """
        touch ${id}.tsv
        """
}

process ANNOTATE_GC {
    label 'gatk_small'

    input:
        tuple val(prefix), path(interval_list)

    output:
        tuple val(prefix), file("${prefix}.annotated.tsv"), emit: annotated_intervals

    script:
        """
        gatk AnnotateIntervals -L $interval_list -R ${params.genome_file} -imr OVERLAPPING_ONLY -O ${prefix}.annotated.tsv
        """

    stub:
        """
        touch ${prefix}.annotated.tsv
        """
}

process CORRECT_GC {

    input:
        tuple val(prefix), file(interval_list), file(annotated)
        tuple val(prefix), val(id), file(tsvs)

    output:
        tuple val(prefix), file("${prefix}.preprocessed.blacklisted.gcfiltered.interval_list"), emit: corrected_intervals
    
    script:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        gatk FilterIntervals -L $interval_list --annotated-intervals $annotated -imr OVERLAPPING_ONLY -O ${prefix}.preprocessed.blacklisted.gcfiltered.interval_list \\
        $tsv_list
        """

    stub:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        touch ${prefix}.preprocessed.blacklisted.gcfiltered.interval_list
        """
}

process REMOVE_ALT_CONTIGS {
    label 'gatk_small'

    input:
        tuple val(prefix), path(interval_list)

    output:
        tuple val(prefix), path("${prefix}.preprocessed.blacklisted.gcfiltered.noalt.interval_list"), emit: noaltcontigs

    script:
        """
        grep -vP "^[A-z,0-9]+_" $interval_list > ${prefix}.preprocessed.blacklisted.gcfiltered.noalt.interval_list
        """

    stub:
        """
        touch ${prefix}.preprocessed.blacklisted.gcfiltered.noalt.interval_list
        """
}

process SCATTER_INTERVALS {
    label 'gatk_small'

    input:
        tuple val(prefix), path(interval_list)
        val(size)

    output:
        tuple val(prefix), path("scatter/*"), emit: scatter

    script:
        count = 0
        """
        mkdir scatter
        gatk IntervalListTools --INPUT $interval_list --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_CONTENT $size --OUTPUT scatter
        ls scatter/*/scattered.interval_list | wc -l
        """

    stub:
        """
        mkdir scatter
        mkdir scatter/temp_001
        mkdir scatter/temp_002
        touch scatter/temp_001/scattered.interval_list
        touch scatter/temp_002/scattered.interval_list
        """
}

process COHORT_PLOIDY {

    input:
        tuple val(prefix), file(interval_list)
        tuple val(prefix), val(id), file(tsvs)

    output:
        tuple val(prefix), path("${prefix}_ploidy-model"), path("${prefix}_ploidy-calls/"), emit: ploidy
    
    script:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        gatk DetermineGermlineContigPloidy \\
        -L $interval_list \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --contig-ploidy-priors ${params.priors} \\
        --output . \\
        --output-prefix ${prefix}_ploidy \\
        $tsv_list
        """

    stub:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        mkdir ${prefix}_ploidy-model ${prefix}_ploidy-calls
        """

}

process COHORT_CALL {

    input:
        tuple val(prefix), val(id), file(tsvs), path(ploidy_model), path(ploidy_calls), file(scatter)

    output:
        tuple val(prefix), path("cohort_calls/cohort_${name}"), emit: calls
    
    script:
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
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output cohort_calls \\
            --output-prefix cohort_${name} \\
        $tsv_list
        touch test
        """

    stub:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        mkdir cohort_calls
        mkdir cohort_calls/cohort_${name}
        """

}

process COHORT_CALL_PANEL {

    input:
        tuple val(prefix), val(id), file(tsvs), path(ploidy_model), path(ploidy_calls), file(intervals)

    output:
        tuple val(prefix), path("cohort_calls/${prefix}*"), emit: calls
    
    script:
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
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output cohort_calls \\
            --output-prefix ${prefix} \\
        $tsv_list
        touch test
        """

    stub:
        tsv_list = tsvs.collect {'-I ' + it}
        tsv_list = tsv_list.join(' ')
        """
        echo $tsv_list
        mkdir cohort_calls
        mkdir cohort_calls/${prefix}
        """

}