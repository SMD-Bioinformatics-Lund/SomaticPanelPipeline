process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    container = "/fs1/resources/containers/multiqc_1_14_0.img"

    input:
        path versions

    output:
        path "software_versions.yml"    , emit: yml
        path "software_versions_mqc.yml", emit: mqc_yml
        path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        template 'dumpsoftwareversions.py'
}
