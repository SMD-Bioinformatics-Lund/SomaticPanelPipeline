process PREPARE_BED {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(vc), path(beds)

    output:
        tuple val(group), val(meta), val(vc), path("out/${task.ext.prefix ?: meta.id}.${task.ext.suffix ?: 'base.bed'}"), emit: base_bed
        path 'versions.yml', emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'base.bed'
        def bed_inputs = beds instanceof List ? beds : [beds]
        """
        mkdir -p out
        python3 ${projectDir}/bin/gens_data.py prepare --beds ${bed_inputs.join(' ')} --out out/${prefix}.${suffix} ${args}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'base.bed'
        """
        mkdir -p out
        touch out/${prefix}.${suffix}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """
}

process EXPAND_ALLELES {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(vc), path(bed)

    output:
        tuple val(group), val(meta), val(vc), path("out/${task.ext.prefix ?: meta.id}.${task.ext.suffix ?: 'expanded.bed'}"), emit: expanded_bed
        path 'versions.yml', emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'expanded.bed'
        """
        mkdir -p out
        python3 ${projectDir}/bin/gens_data.py expand --bed ${bed} --out out/${prefix}.${suffix} ${args}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'expanded.bed'
        """
        mkdir -p out
        touch out/${prefix}.${suffix}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """
}

process WRITE_COMMANDS {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(baf), path(baf_index), path(cov), path(cov_index)

    output:
        tuple val(group), val(meta), path("${task.ext.prefix ?: meta.id}.gens"), emit: dbload
        tuple val(group), val(meta), path("${task.ext.prefix ?: meta.id}.gens_v4_somatic"), emit: gens_v4
        path 'versions.yml', emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def process_group = meta.paired ? "${group}p" : "${group}"
        """
        python3 ${projectDir}/bin/gens_data.py commands \\
            --sample-id ${meta.id} \\
            --case-id ${process_group} \\
            --sample-type ${meta.type} \\
            --access-dir ${params.gens_accessdir} \\
            --baf ${baf.name} \\
            --coverage ${cov.name} \\
            --prefix ${prefix} \\
            ${args}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.gens ${prefix}.gens_v4_somatic
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """
}
