process SEQTK {
	cpus 10
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
		tuple val(group), val(meta), file("*R1*fastq.gz"), file("*R2*fastq.gz"), emit: fastq_sub

    script:
		"""
		seqtk sample -s 1234 $r1 ${meta.sub}  | gzip --no-name > ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R1_subsample${meta.sub}.fastq.gz &
		seqtk sample -s 1234 $r2 ${meta.sub}  | gzip --no-name > ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R2_subsample${meta.sub}.fastq.gz
		"""
    stub:
		"""
		touch ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R1_subsample${meta.sub}.fastq.gz
		touch ${meta.clarity_sample_id}_${meta.clarity_pool_id}_R2_subsample${meta.sub}.fastq.gz
		"""

}