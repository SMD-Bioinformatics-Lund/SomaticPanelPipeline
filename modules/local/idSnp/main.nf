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
        bcftools mpileup $args $bam | bcftools call $args2 > ${prefix}.vcf
        
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
        tuple val(group), val(meta), file(vcfs)

    output:
        tuple val(group), file("*.csv"),            emit: idsnp_checked
        path "versions.yml", optional: true         ,emit:   versions

    when:
        task.ext.when == null || task.ext.when

    script:
        
        def args    = task.ext.args  ?: ""
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
	    normalvcf   = vcfs[normal_idx]
	    tumorvcf    = vcfs[tumor_idx]

        if(meta.id.size() == 2) {
            """
            idsnp_controller-myeolid_specific.pl \\
                --vcf_sample $tumorvcf  \\
                --vcf_control $normalvcf \\
                --sample  $tumor_id \\
                --control $normal_id \\
                $args
    
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
            """
            echo "Not applicable" > s${tumor_id}.csv
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
        
        if(meta.id.size() == 2) {
            """ 
            touch s${tumor_id}_c${normal_id}.csv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
        
            """
            echo "Not applicable" > s${tumor_id}.csv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process PROVIDER {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    
    output:
        tuple val(group), val(meta), file("*.genotypes"),    emit: genotype_checked
        path "versions.yml",                                 emit: versions  // Emit version information in YAML format
    
    when:
        task.ext.when == null || task.ext.when
        
    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""

        """
        provider.pl  --bam $bam  $args  --out $prefix 
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
        
    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"

        """
        touch ${prefix}.genotypes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}