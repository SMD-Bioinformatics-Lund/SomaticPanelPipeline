process CNVKIT2ASCAT {
	cpus 1
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(id), val(type), file(baflogr)

	output:
		tuple val(group), val(id), val(type), file("tumor_baf.txt"), file("tumor_log2r.txt"), emit: ascat_input

	"""
	cnvkit2ascat.pl $baflogr $id
	"""
}



process ASCAT_252 {
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/", mode: 'copy', overwrite: true, pattern: '*.txt'
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	container = '/fs1/resources/containers/ascat2.5.2_patched.sif'

	input:
		tuple val(group), val(id), val(type), file(baf), file(logr2)

	output:
		tuple val(group), val(id), val(type), file("HRD/${group}_HRD_Telli.txt"), emit: ascat252_HRD_score
		tuple val(group), val(id), val(type), file("*.png")

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/run_ascat.R 
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/GenomicScars.R $params.scarfunctions $params.chrominfo_solid $params.chrominfo_general $params.fdata $group
	"""

}


process ASCAT_30 {
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/", mode: 'copy', overwrite: true, pattern: '*.txt'
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	container = '/fs1/resources/containers/ascat3.0_patched.sif'

	input:
		tuple val(group), val(id), val(type), file(baf), file(logr2)

	output:
		tuple val(group), val(id), val(type), file("HRD/${group}_HRD_Telli.txt"), emit: ascat30_HRD_score
		tuple val(group), val(id), val(type), file("*.png")

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/run_ascat.R
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/GenomicScars.R $params.scarfunctions $params.chrominfo_solid $params.chrominfo_general $params.fdata $group
	"""

}

