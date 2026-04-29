process VALIDATE_COYOTE_SNV {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        path(assay_config)
        path(known_variants)
        
    output:
        tuple val(group), val(meta), file("*results.txt"),      emit: results
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix   ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        simulate_coyote_default_filters.py --vcf $vcf --known $known_variants --config $assay_config > ${prefix}.results.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
    stub:
        def prefix  = task.ext.prefix ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        touch ${prefix}.results.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
        
}
