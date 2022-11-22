process QC_TO_CDM {
	cpus 1
	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	time '10m'
	memory '50 MB'

	input:
		tuple val(id), val(type), file(qc), val(clarity_sample_id), val(sequencing_run), val(r1), val(r2)

	output:
		tuple val(id), file("${id}.cdm"), emit: cdm_done

	when:
		!params.noupload

	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")

	"""
	echo "--run-folder $rundir --sample-id $id --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/${id}_${type}.QC --lims-id $clarity_sample_id" > ${id}.cdm
	"""
}

process LOWCOV {
	cpus 1
	memory '5 GB'
	publishDir "${params.outdir}/${params.subdir}/QC", mode: 'copy', overwrite: 'true'
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(id), val(type), file(bam), file(bai), file(dedup) //from bam_lowcov


	output:
		tuple val(group), val(id), val(type), file("${id}.lowcov.bed"), emit: lowcov_regions

	"""
    source activate sambamba
	panel_depth.pl $bam ${params.regions_proteincoding} > lowcov.bed
	overlapping_genes.pl lowcov.bed ${params.gene_regions} > ${id}.lowcov.bed
	"""
}