process CSV_CHECK {
    label 'process_single'

    input:
        path samplesheet

    output:
        path "*.checked.csv", emit: csv

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${samplesheet.baseName}"
        """
        check_samplesheet.py -c ${samplesheet} -o samplecheck.txt --cmd-assays-json ${projectDir}/resources/cmd_assays.json

        if [[ -e "samplecheck.txt" ]]; then
            cp ${samplesheet} "${prefix}.checked.csv"
        else
            echo "samplecheck.txt does not exist"
        fi
        """
}

process UNSPRING {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(spring)

    output:
        tuple val(group), val(meta), path("*.unspring_R1.fastq.gz"), path("*.unspring_R2.fastq.gz"), emit: fastq

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        """
        spring -d \\
            -i ${spring} \\
            -o ${prefix}.unspring_R1.fastq.gz ${prefix}.unspring_R2.fastq.gz \\
            -g \\
            -t ${task.cpus} $args
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.type}"
        """
        touch ${prefix}.unspring_R1.fastq.gz
        touch ${prefix}.unspring_R2.fastq.gz
        """
}
