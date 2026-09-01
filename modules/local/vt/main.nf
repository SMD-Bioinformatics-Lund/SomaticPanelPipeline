process VT_DECOMPOSE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(vcf), path(tbi)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf.gz"), emit: decomposed_vcf
        path "versions.yml",                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "decomposed.vcf.gz"
        """
        mkdir -p out
        vt decompose $args $vcf -o out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "decomposed.vcf.gz"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
        END_VERSIONS
        """
}

process VT_NORMALIZE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf.gz"), emit: normalized_vcf
        path "versions.yml",                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "normalized.vcf.gz"
        """
        mkdir -p out
        vt normalize $vcf $args -o out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "normalized.vcf.gz"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """
}

process VT_UNIQ {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf.gz"), emit: unique_vcf
        path "versions.yml",                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "vcf.gz"
        """
        mkdir -p out
        vt uniq $args $vcf -o out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-uniq: \$(echo \$(vt uniq 2>&1) | sed 's/.*uniq v//; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}_${vc}"
        def suffix = task.ext.suffix ?: "vcf.gz"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-uniq: \$(echo \$(vt uniq 2>&1) | sed 's/.*uniq v//; s/ .*//')
        END_VERSIONS
        """
}
