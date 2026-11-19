process ANNOTATE_VEP {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), path(vcf)

    output:
        tuple val(group), val(meta), path("*.vep.vcf"), emit: vcf_vep
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${vcf.baseName}"

        """
        vep -i ${vcf} -o ${prefix}".vep.vcf" --fork ${task.cpus} $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${vcf.baseName}"

        """
        touch ${prefix}".vep.vcf"
        echo $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        END_VERSIONS
        """
}
