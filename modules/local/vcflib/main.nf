process FILTER {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), file("out/*.vcf"), emit: filtered_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "filtered.vcf"

        """
        mkdir -p out
        vcffilter $args $vcf > out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$(vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "filtered.vcf"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$(vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g')
        END_VERSIONS
        """
}

process GLXGT {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), file("out/*.vcf"),    emit: fixed_genotypes_vcf
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "glxgt.vcf"

        """
        mkdir -p out
        vcfglxgt < $vcf > out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfglxgt: \$(echo \$(vcfglxgt -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "glxgt.vcf"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfglxgt: \$(echo \$(vcfglxgt -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g')
        END_VERSIONS
        """
}
