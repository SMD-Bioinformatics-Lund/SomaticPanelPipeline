process BCFTOOLS_NORM {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(vcf), path(tbi)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf.gz"), path("out/*.vcf.gz.tbi"), emit: normalized_vcfs
        path "versions.yml",                                                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "norm.vcf.gz"
        def output_vcf = "${prefix}.${suffix}"

        """
        mkdir -p out
        bcftools norm ${vcf} ${args} -o out/${output_vcf}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "norm.vcf.gz"
        def output_vcf = "${prefix}.${suffix}"
        """
        mkdir -p out
        touch out/${output_vcf}
        touch out/${output_vcf}.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """
}

process BCFTOOLS_MPILEUP_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(bam), path(bai)

    output:
        tuple val(group), val(meta), path("*.raw.vcf"), emit: raw_vcf
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        """
        bcftools mpileup $args $bam | bcftools call $args2 > ${prefix}.raw.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.raw.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """
}

process BCFTOOLS_QUERY_IDSNP {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(vcf)

    output:
        tuple val(group), val(meta), path(vcf), path("*.genotypes"), emit: genotypes
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        bcftools query -f '%ID\\t[%GT]\\n' ${vcf} > ${prefix}.genotypes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.genotypes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """
}

process ANNOTATE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf"), emit: annotated_vcf
        path "versions.yml",                                     emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "annotated.vcf"

        """
        mkdir -p out
        bcftools annotate $args $vcf -o out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        def suffix = task.ext.suffix ?: "annotated.vcf"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | sed -n '1{ s/^bcftools //; s/ .*//; p }')
        END_VERSIONS
        """
}
