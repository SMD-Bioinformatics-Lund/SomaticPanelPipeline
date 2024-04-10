process ALLELE_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.vcf"), emit:   sample_id_vcf
        path "versions.yml",                        emit:   versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        """
        bcftools mpileup \\
                 $args $bam \\ | 
        bcftools call \\
                 $arg2 > ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """
}

process SNP_CHECK {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(tumorvcf), file(normalvcf)

    output:
        tuple val(group), val(meta), file("*.csv"), emit: idsnp_checked

    when:
        task.ext.when == null || task.ext.when

    script:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]

        """
        idsnp_controller-myeolid_specific.pl \\
            --vcf_sample $tumorvcf  \\
            --vcf_control $normalvcf \\
            --sample  $tumor_id \\
            --control $normal_id 
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
        """ 
        touch s$tumor_id_c$normal.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}