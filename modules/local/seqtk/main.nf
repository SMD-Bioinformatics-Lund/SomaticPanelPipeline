process SEQTK {
    cpus 2
    memory '128 GB'
    time '2h'
    errorStrategy 'retry'
    maxErrors 5
    tag "${meta.id}"
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    container = "/fs1/resources/containers/seqtk_1.3.sif"

    input:
        tuple val(group), val(meta), file(r1), file(r2)

    output:
        tuple val(group), val(meta), file("${meta.clarity_sample_id}_${meta.clarity_pool_id}_R1_subsample${meta.sub}.fastq.gz"), file("${meta.clarity_sample_id}_${meta.clarity_pool_id}_R2_subsample${meta.sub}.fastq.gz"),                    emit: fastq_sub
        path "versions.yml",           emit: versions

    script:
        """
        seqtk sample -s 1234 $r1 ${meta.sub} | gzip --no-name > ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R1_subsample${meta.sub}.fastq.gz &
        seqtk sample -s 1234 $r2 ${meta.sub} | gzip --no-name > ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R2_subsample${meta.sub}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R1_subsample${meta.sub}.fastq.gz
        touch ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R2_subsample${meta.sub}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/.*Version: //; s/ .*//')
        END_VERSIONS
        """
}