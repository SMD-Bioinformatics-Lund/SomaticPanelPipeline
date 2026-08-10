process CREATE_SNVPON {
    label "process_alot"
    tag "$vc"

    input:
        tuple val(group), val(vc), path(vcfs)

    output:
        tuple val(group), val(vc), path("*_${vc}_PON.snv"), emit: SNV_PON
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        create_snv_pon.py --vcf-mask "*.vcf.gz" --out ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        echo $vcfs
        touch ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}
