process TABIX_BGZIPTABIX {
    tag "$group"
    label 'process_single'

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("*.gz"), path("*.{tbi,csi}"),    emit: gz_index
        path "versions.yml",                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"
        """
        bgzip --threads ${task.cpus} -c $args $vcf > ${prefix}.${vcf.getExtension()}.gz
        tabix $args2 ${prefix}.${vcf.getExtension()}.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --version 2>&1 | sed -n '1{ s/^bgzip (htslib) //; s/ .*//; p }')
            tabix: \$(tabix --version 2>&1 | sed -n '1{ s/^tabix (htslib) //; s/ .*//; p }')
        END_VERSIONS

        """

    stub:
        def args2 = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"
        def index = args2.contains("-C ") || args2.contains("--csi") ? "csi" : "tbi"
        """
        echo "" | gzip > ${prefix}.${vcf.getExtension()}.gz
        touch ${prefix}.${vcf.getExtension()}.gz.${index}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --version 2>&1 | sed -n '1{ s/^bgzip (htslib) //; s/ .*//; p }')
            tabix: \$(tabix --version 2>&1 | sed -n '1{ s/^tabix (htslib) //; s/ .*//; p }')
        END_VERSIONS
    """
}

process TABIX_INDEX {
    tag "$group"
    label 'process_single'

    input:
        tuple val(group), val(meta), val(vc), path(input_file)

    output:
        tuple val(group), val(meta), val(vc), path("out/${input_file.name}"), path("out/${input_file.name}.{tbi,csi}"), emit: indexed_file
        path "versions.yml",                                                              emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def index = args.contains("-C ") || args.contains("--csi") ? "csi" : "tbi"
        """
        mkdir -p out
        tabix $args $input_file
        ln -sf ../$input_file out/${input_file.name}
        ln -sf ../${input_file}.${index} out/${input_file.name}.${index}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: \$(tabix --version 2>&1 | sed -n '1{ s/^tabix (htslib) //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def index = args.contains("-C ") || args.contains("--csi") ? "csi" : "tbi"
        """
        mkdir -p out
        touch out/${input_file.name}
        touch out/${input_file.name}.${index}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: \$(tabix --version 2>&1 | sed -n '1{ s/^tabix (htslib) //; s/ .*//; p }')
        END_VERSIONS
        """
}

process BGZIP {
    label 'process_single'
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), path(input_file)

    output:
        tuple val(group), val(meta), val(vc), path("out/${task.ext.prefix ?: input_file.baseName}.${task.ext.suffix ?: (input_file.extension + '.gz')}"), emit: compressed_file
        path 'versions.yml', emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${input_file.baseName}"
        def suffix = task.ext.suffix ?: "${input_file.extension}.gz"
        """
        mkdir -p out
        bgzip -c ${args} ${input_file} > out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --version 2>&1 | sed -n '1{ s/^bgzip (htslib) //; s/ .*//; p }')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${input_file.baseName}"
        def suffix = task.ext.suffix ?: "${input_file.extension}.gz"
        """
        mkdir -p out
        touch out/${prefix}.${suffix}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --version 2>&1 | sed -n '1{ s/^bgzip (htslib) //; s/ .*//; p }')
        END_VERSIONS
        """
}
