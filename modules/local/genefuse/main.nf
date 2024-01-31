process GENEFUSE {
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("*.genefuse.json"),   emit: genefuse_json
        tuple val(group), val(meta), file("*.html"),            emit: genefuse_html
        path "versions.yml",                                    emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        genefuse -t $task.cpus \\
            $args \\
            -1 $r1 \\
            -2 $r2 \\
            -h ${prefix}.html \\
            -j ${prefix}.genefuse.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            genefuse: 0.8.0
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.genefuse.json 
        touch ${prefix}.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            genefuse: 0.8.0
        END_VERSIONS
        """
}