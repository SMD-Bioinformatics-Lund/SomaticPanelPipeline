process CREATE_SNVPON {
    label "process_alot"
    tag "$vc"

    input:
        tuple val(group), val(meta), val(vc), path(vcfs), path(tbis)

    output:
        tuple val(group), val(meta), val(vc), path("*_${vc}_PON.snv"), emit: SNV_PON
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        create_snv_pon.pl "*.vcf.gz" > ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${params.assay}"
        """
        echo $vcfs
        touch ${prefix}_${vc}_PON.snv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}
