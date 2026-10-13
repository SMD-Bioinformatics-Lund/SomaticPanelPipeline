process VCFANNO {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), path(vcf)

    output:
        tuple val(group), val(meta), path("*.agg.enigma.vcf"), emit: vcf_enigma
        path "versions.yml",                                   emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        vcfanno $args $vcf > ${prefix}.agg.enigma.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfanno: \$( echo \$(vcfanno 2>&1) | sed 's/.*version //' | sed 's/ \\[.*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.agg.enigma.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcfanno: \$( echo \$(vcfanno 2>&1) | sed 's/.*version //' | sed 's/ \\[.*//')
        END_VERSIONS
        """
}
