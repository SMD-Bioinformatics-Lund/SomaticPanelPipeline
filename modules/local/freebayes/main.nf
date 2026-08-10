process FREEBAYES {
    label "process_low"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        each file(bed)

    output:
        tuple val(group), val(meta), val("freebayes"), file("freebayes_${bed}.vcf.raw"),    emit: vcfparts_freebayes
        path "versions.yml",                                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''

        if( meta.id.size() >= 2 ) {

            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            freebayes -t $bed \\
            $args \\
            -F ${params.fb_var_freq_cutoff_p} \\
            ${bams[tumor_idx]} \\
            ${bams[normal_idx]} > freebayes_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            freebayes -t $bed \\
            $args \\
            -F ${params.fb_var_freq_cutoff_up} \\
            $bams > freebayes_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            END_VERSIONS
            """
        }

    stub:
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            echo tumor:${bams[tumor_idx]} ${meta.id[tumor_idx]} normal:${bams[normal_idx]} ${meta.id[normal_idx]}
            touch freebayes_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            END_VERSIONS
            """
        }
        else {
            """
            echo tumor:$bams
            touch freebayes_${bed}.vcf.raw

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            END_VERSIONS
            """
        }
}
