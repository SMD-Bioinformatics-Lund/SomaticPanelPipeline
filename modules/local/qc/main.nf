process QC_TO_CDM {
	cpus 1
	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: 'true'
	tag "${meta.id}"
	time '10m'
	memory '50 MB'

	input:
		tuple val(group), val(meta), file(qc)

	output:
		tuple val(group), val(meta), file("${meta.id}.cdm"), emit: cdm_done

	when:
		!params.noupload

	script:
	"""
	echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${meta.id}.cdm
	"""
	stub:
	"""
	echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${meta.id}.cdm
	"""
}

process LOWCOV {
	cpus 1
	memory '5 GB'
	publishDir "${params.outdir}/${params.subdir}/QC", mode: 'copy', overwrite: 'true'
	time '1h'
	tag "${meta.id}"

	input:
		tuple val(group), val(meta), file(bam), file(bai), file(dedup) //from bam_lowcov

	output:
		tuple val(group), val(meta.type), file("${meta.id}.lowcov.bed"), emit: lowcov_regions

	script:
	"""
    source activate sambamba
	panel_depth.pl $bam ${params.regions_proteincoding} > lowcov.bed
	overlapping_genes.pl lowcov.bed ${params.gene_regions} > ${meta.id}.lowcov.bed
	"""
	stub:
	"""
	touch ${meta.id}.lowcov.bed
	"""
}