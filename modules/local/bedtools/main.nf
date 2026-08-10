process BEDTOOLS_INTERSECT {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(vc), file(input_file)
        path(intersect)

    output:
        tuple val(group), val(meta), file("${out_file}"),          emit: intersected
        tuple val(group), val(meta), val(vc), file("${out_file}"), emit: vcf_intersected
        path "versions.yml",                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args ?: ''
        def prefix      = task.ext.prefix ?: "${meta.id}"
        def suffix      = task.ext.suffix ?: 'bed'
        out_file        = "${prefix}.${suffix}"
        """
        bedtools intersect -a $input_file -b $intersect $args > ${out_file}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: "${meta.id}"
        def suffix      = task.ext.suffix ?: 'bed'
        out_file        = "${prefix}.${suffix}"
        """
        echo $input_file
        touch ${out_file}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
}
