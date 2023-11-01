process MSISENSOR {
    label 'process_low'
    label 'error_ignore'
    tag "$group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
        
    output:
        tuple val(group), file("*.msi_single"),                 emit: msi_score
        tuple val(group), file("*.msi_paired"), optional:true,  emit: msi_score_paired
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''
        def prefix  = task.ext.prefix   ?: "${group}"

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            msisensor-pro msi \\
                $args \\
                -n ${bams[normal_idx]} \\
                -t ${bams[tumor_idx]} \\
                -o ${prefix}.msi_paired
            
            msisensor-pro pro \\
                $args2 \\
                -t ${bams[tumor_idx]} \\
                -o ${prefix}.msi_single

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                msisensor-pro: \$(echo \$(msisensor-pro --version 2>&1) | grep 'Version: v' | sed 's/Version: v//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            msisensor-pro pro \\
                $args2 \\
                -t $bams \\
                -o ${prefix}.msi_single

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                msisensor-pro: \$(echo \$(msisensor-pro --version 2>&1) | grep 'Version: v' | sed 's/Version: v//')
            END_VERSIONS
            """
        }

    stub:
        def prefix  = task.ext.prefix   ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            touch ${prefix}.msi_single

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                msisensor-pro: \$(echo \$(msisensor-pro --version 2>&1) | grep 'Version: v' | sed 's/Version: v//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            touch ${prefix}.msi_single

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                msisensor-pro: \$(echo \$(msisensor-pro --version 2>&1) | grep 'Version: v' | sed 's/Version: v//')
            END_VERSIONS
            """
        }
}