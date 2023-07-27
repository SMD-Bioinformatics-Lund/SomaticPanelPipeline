process MELT {
    label "process_low"
	tag "${meta.id}"

	when:
		params.melt

	input:
		tuple val(group), val(meta), file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

	output:
		tuple val(group), val("melt"), file("${meta.id}.melt.merged.vcf"), emit: melt_vcf

    script:
        """
        java -jar /opt/MELTv2.2.2/MELT.jar Single \\
            -bamfile $bam \\
            -r 150 \\
            -h $params.genome_file \\
            -n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
            -z 50000 \\
            -d 50 -t /opt/mei_list \\
            -w . \\
            -c $MEAN_DEPTH \\
            -cov $COV_DEV \\
            -e $INS_SIZE
        merge_melt.pl $params.meltheader ${meta.id}
        """
    stub:
        """
        touch ${meta.id}.melt.merged.vcf
        """
}