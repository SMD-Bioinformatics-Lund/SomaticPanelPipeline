process GENS_YAML {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta_list), file(baf_files), file(cov_files), file(loh_cat)

    output:
        tuple val(group), file("*.gens_somatic.yaml"), emit: gens_yaml
        path "versions.yml",                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix        = task.ext.prefix ?: "${group}"
        def process_group = meta_list.any { it.paired } ? group + 'p' : group
        def sample_ids    = meta_list.collect { it.id }
        def sample_types  = meta_list.collect { it.type }
        def sexes         = meta_list.collect { it.sex ?: '.' }
        def baf_names    = [baf_files].flatten().collect { it.name }
        def cov_names    = [cov_files].flatten().collect { it.name }

        """
        gens_yaml.py \\
            --case_id ${process_group} \\
            --genome_build 38 \\
            --gens_accessdir ${params.gens_accessdir} \\
            --sample_ids ${sample_ids.join(' ')} \\
            --sample_types ${sample_types.join(' ')} \\
            --sample_sexes ${sexes.join(' ')} \\
            --baf_filenames ${baf_names.join(' ')} \\
            --cov_filenamess ${cov_names.join(' ')} \\
            --loh_cat_filename ${loh_cat} \\
            --output ${prefix}.gens_somatic.yaml

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.gens_somatic.yaml

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}



