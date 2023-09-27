process FASTP {
    cpus 10
    memory '128 GB'
    time '2h'
    errorStrategy 'retry'
    maxErrors 5
    tag "${meta.id}"
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    container = "/fs1/resources/containers/fastp_0.23.sif"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("${meta.id}.${meta.type}.cleaned_R1.fastq.gz"), file("${meta.id}.${meta.type}.cleaned_R2.fastq.gz"),	emit: fastq_trimmed
        tuple val(group), file("${meta.id}.${meta.type}.cleaned.fastp.json"), file("${meta.id}.${meta.type}.cleaned.html"), 					emit: fastq_stats
        path "versions.yml", 																													emit: versions

    script:
        """
        fastp -w ${task.cpus} \\
        --average_qual 20 --length_required 50 --qualified_quality_phred 20 \\
        -i $r1 -I $r2 -o ${meta.id}.${meta.type}.cleaned_R1.fastq.gz -O ${meta.id}.${meta.type}.cleaned_R2.fastq.gz \\
        -h ${meta.id}.${meta.type}.cleaned.html -j ${meta.id}.${meta.type}.cleaned.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.${meta.type}.cleaned_R1.fastq.gz
        touch ${meta.id}.${meta.type}.cleaned_R2.fastq.gz
        touch ${meta.id}.${meta.type}.cleaned.html
        touch ${meta.id}.${meta.type}.cleaned.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
}