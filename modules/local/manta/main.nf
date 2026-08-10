process MANTA {
    label 'process_medium'
    label 'stage'
    label 'scratch'
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(bam), file(bai)
        val(reference)
        val(type)

    output:
        tuple val(group), val(meta), file("*_manta.*.vcf"), emit: manta_vcf
        path "versions.yml",                                emit: versions
        

    when:
        task.ext.when == null || task.ext.when
    
    script:
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal      = normal_idx >= 0 ? bam[normal_idx] : null
        normal_id   = normal_idx >= 0 ? meta.id[normal_idx] : ''
        tumor       = bam[tumor_idx]
        tumor_id    = meta.id[tumor_idx]
        prefix      = task.ext.prefix  ?: tumor_id
        prefix2     = task.ext.prefix2 ?: normal_id

        if(meta.id.size() == 2) { 
            """
            set +eu
            source activate py2
            set -eu
            configManta.py \\
                --tumorBam $tumor \\
                --normalBam $normal \\
                --callRegions $reference \\
                $args \\
                --runDir .
            python runWorkflow.py $args2
            mv results/variants/somaticSV.vcf.gz ${prefix}_manta.${type}.vcf.gz
            mv results/variants/diploidSV.vcf.gz ${prefix2}_manta.${type}.vcf.gz
            gunzip ${prefix}_manta.${type}.vcf.gz
            gunzip ${prefix2}_manta.${type}.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else {
            """
            set +eu
            source activate py2
            set -eu
            configManta.py \\
                --tumorBam $tumor \\
                --callRegions $reference \\
                $args \\
                --runDir .
            python runWorkflow.py $args2
            mv results/variants/tumorSV.vcf.gz ${prefix}_manta.${type}.vcf.gz
            gunzip ${prefix}_manta.${type}.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal      = normal_idx >= 0 ? bam[normal_idx] : null
        normal_id   = normal_idx >= 0 ? meta.id[normal_idx] : ''
        tumor       = bam[tumor_idx]
        tumor_id    = meta.id[tumor_idx]
        prefix      = task.ext.prefix  ?: tumor_id
        prefix2     = task.ext.prefix2 ?: normal_id
        if(meta.id.size() == 2) {
            """
            set +eu
            source activate py2
            set -eu
            touch ${prefix}_manta.${type}.vcf 
            touch ${prefix2}_manta.${type}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else {
            """
            set +eu
            source activate py2
            set -eu
            touch ${prefix}_manta.${type}.vcf 
            touch ${prefix}_manta_filtered.${type}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                manta: \$(configManta.py --version)
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

}
