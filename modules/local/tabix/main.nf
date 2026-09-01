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
        tuple val(group), val(meta), val(vc), path(vcf)

    output:
        tuple val(group), val(meta), val(vc), path("out/*.vcf.gz"), path("out/*.{tbi,csi}"), emit: indexed_vcf
        path "versions.yml",                                                              emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def index = args.contains("-C ") || args.contains("--csi") ? "csi" : "tbi"
        """
        mkdir -p out
        tabix $args $vcf
        ln -sf ../$vcf out/${vcf.name}
        ln -sf ../${vcf}.${index} out/${vcf.name}.${index}

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
        touch out/${vcf.name}
        touch out/${vcf.name}.${index}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: \$(tabix --version 2>&1 | sed -n '1{ s/^tabix (htslib) //; s/ .*//; p }')
        END_VERSIONS
        """
}
