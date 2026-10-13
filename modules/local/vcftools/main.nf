process VCF_CONCAT {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcfs)

    output:
        tuple val(group), val(meta), val(vc), file("*.concat.vcf"),     emit: concatenated_vcf
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${vc}"
        """
        vcf-concat $vcfs | vcf-sort -c > ${prefix}.concat.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}_${vc}_bcftools.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
        END_VERSIONS
        """
}

process VCF_SORT {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf"), emit: sorted_vcf
        path "versions.yml",                                     emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: (vc ? "${vc}" : "${group}")
        def suffix = task.ext.suffix ?: "sorted.vcf"
        """
        mkdir -p out
        vcf-sort -c $vcf > out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: (vc ? "${vc}" : "${group}")
        def suffix = task.ext.suffix ?: "sorted.vcf"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
        END_VERSIONS
        """
}
