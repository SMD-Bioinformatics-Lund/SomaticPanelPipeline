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