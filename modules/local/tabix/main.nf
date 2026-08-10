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
            bgzip: \$(echo \$(bgzip --version) | sed '1!d;s/.* //')
            tabix: \$(echo \$(tabix -h 2>&1) | grep -oP 'Version:\\s*\\K[^\\s]+')
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
            bgzip: \$(echo \$(bgzip --version) | sed '1!d;s/.* //')
            tabix: \$(echo \$(tabix -h 2>&1) | grep -oP 'Version:\\s*\\K[^\\s]+')
        END_VERSIONS
    """
}
