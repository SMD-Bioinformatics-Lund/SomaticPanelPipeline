process MELT {
    label "process_medium", "error_retry"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

    output:
        tuple val(group), val(meta), val("melt"), file("*.melt.merged.vcf"),    emit: melt_vcf
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        java -jar /opt/MELTv2.2.2/MELT.jar Single \\
            -bamfile $bam \\
            $args \\
            -w . \\
            -c $MEAN_DEPTH \\
            -cov $COV_DEV \\
            -e $INS_SIZE

        merge_melt.pl ${params.meltheader} ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            melt: \$(java -jar /opt/MELTv2.2.2/MELT.jar 2>&1 | grep ^MELTv | sed 's/MELTv//' | cut -d' ' -f1)
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.melt.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            melt: \$(java -jar /opt/MELTv2.2.2/MELT.jar 2>&1 | grep ^MELTv | sed 's/MELTv//' | cut -d' ' -f1)
        END_VERSIONS
        """
}