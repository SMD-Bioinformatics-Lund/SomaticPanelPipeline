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
	publishDir "${params.outdir}/${params.subdir}/ASCAT2.5.2/", mode: 'copy', overwrite: true, pattern: '*.txt'
	publishDir "${params.outdir}/${params.subdir}/ASCAT2.5.2/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	errorStrategy 'ignore'
	container = '/fs1/resources/containers/ascat2.5.2_patched.sif'

	input:
		tuple val(group), val(id), val(type), file(baf), file(logr2)

	output:
		tuple val(group), val(id), val(type), file("${group}_HRD_Telli.txt"), emit: ascat252_HRD_score
		tuple val(group), val(id), val(type), file("*.png"), emit: ascat252_plots
		tuple val(group), val(id), val(type), file("${group}.results.txt"), emit: ascat252_purplo

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/run_ascat.R 
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/GenomicScars.R $params.scarfunctions $params.chrominfo_solid $params.chrominfo_general $params.fdata $group
	cp HRD/${group}_HRD_Telli.txt .
	mv ASCAT_summary_results_gamma0.9_penalty20.txt ${group}.results.txt
	"""

}


process ASCAT_30 {
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/", mode: 'copy', overwrite: true, pattern: '*.txt'
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	errorStrategy 'ignore'
	container = '/fs1/resources/containers/ascat3.0_patched.sif'

	input:
		tuple val(group), val(id), val(type), file(baf), file(logr2)

	output:
		tuple val(group), val(id), val(type), file("${id}_HRD_Telli.txt"), emit: ascat30_HRD_score
		tuple val(group), val(id), val(type), file("*.png"), emit: ascat30_plots
		tuple val(group), val(id), val(type), file("${id}.results.txt"), emit: ascat30_purplo
		tuple val(group), val(id), val(type), file("${id}.LogR.PCFed.txt"), file("${id}.BAF.PCFed.txt"), emit: baflogr
		tuple val(group), val(id), file("${id}.ploidy"), emit: ploidy

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/run_ascat_3.0.R
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/workflows/bin/GenomicScars.R $params.scarfunctions $params.chrominfo_solid $params.chrominfo_general $params.fdata $id
	cp HRD/${id}_HRD_Telli.txt .
	mv ASCAT_summary_results_gamma0.9_penalty20.txt ${id}.results.txt
	grep -v ^Sample ${id}.results.txt | cut -f 4 > ${id}.ploidy
	"""

}

