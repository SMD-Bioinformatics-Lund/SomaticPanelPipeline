process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    input:
        path versions
        tuple val(group), val(meta)

    output:
        path "*.software_versions.yml"    , emit: yml
        path "*.software_versions_mqc.yml", emit: mqc_yml
        path "*.versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        prefix      = task.ext.prefix ?: "${meta.group}"
        template 'dumpsoftwareversions.py'

}
