process INTERSECT {
    label "process_low"
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf)
        val(intersect)
        
    output:
        tuple val(group), val(meta), file("${prefix}.intersected.${suffix}"), emit: intersected
        path "versions.yml",                                      emit: versions

    script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        suffix = task.ext.suffix ?: ''
        """
        bedtools intersect $args -a $vcf -b $intersect > ${prefix}.intersected.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools | grep Version | sed -r "s/Version:\s+//")
        END_VERSIONS
        """
    stub:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        suffix = task.ext.suffix ?: ''
        """
        echo $intersect $args > ${prefix}.intersected.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools | grep Version | sed -r "s/Version:\s+//")
        END_VERSIONS
        """
}   