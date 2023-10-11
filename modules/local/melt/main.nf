process MELT {
    label "process_low"
    tag "${meta.id}"

    when:
        params.melt

    input:
        tuple val(group), val(meta), file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

    output:
        tuple val(group), val(meta), val("melt"), file("${meta.id}.melt.merged.vcf"),  emit: melt_vcf
        path "versions.yml",                                                emit: versions

    script:
        """
        java -jar /opt/MELTv2.2.2/MELT.jar Single \\
            -bamfile $bam \\
            -r 150 \\
            -h $params.genome_file \\
            -n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
            -z 500000 \\
            -d 50 -t /opt/mei_list \\
            -w . \\
            -c $MEAN_DEPTH \\
            -cov $COV_DEV \\
            -e $INS_SIZE
        merge_melt.pl $params.meltheader ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            melt: \$(java -jar /opt/MELTv2.2.2/MELT.jar 2>&1 | grep ^MELTv | sed 's/MELTv//' | cut -d' ' -f1)
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.melt.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            melt: \$(java -jar /opt/MELTv2.2.2/MELT.jar 2>&1 | grep ^MELTv | sed 's/MELTv//' | cut -d' ' -f1)
        END_VERSIONS
        """
}