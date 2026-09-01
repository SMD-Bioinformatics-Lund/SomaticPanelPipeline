process DEPTH {
    label 'process_medium'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.panel.depth.tsv"), emit: depth_tsv
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ""
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        sambamba depth base $bam $args > ${prefix}.panel.depth.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sambamba: \$(echo \$(sambamba --version 2>&1) | sed 's/sambamba //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat <<-EOF > ${prefix}.panel.depth.tsv
        REF\tPOS\tCOV
        chr1\t100\t600
        chr1\t101\t400
        chr1\t102\t300
        EOF

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sambamba: \$(echo \$(sambamba --version 2>&1) | sed 's/sambamba //; s/ .*//')
        END_VERSIONS
        """
}