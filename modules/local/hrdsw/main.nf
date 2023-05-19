process CNVKIT2OVAHRDSCAR {
	cpus 1
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(id), val(type), file(segments)

	output:
		tuple val(group), val(id), val(type), val(caller), file("${id}.cnvkit.ovaHRDscar.txt"), emit: ovaHRDscar_segments

	script:
		caller = "cnvkit"
		if (segments =~ /purity/) {
			caller = "cnvkitpurity"
		}

	"""
	cnvkit2HRD.pl $segments $id > ${id}.cnvkit.ovaHRDscar.txt
	"""
}

process CNVKIT2SCARHRD {
 	cpus 1
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(meta), val(part), file(segments)

	output:
		tuple val(group), val(meta), val(caller), file("${meta.id}.cnvkit.scarHRD.txt"), emit: scarHRD_segments

	script:
		ploidyv = "NA"
		caller = "cnvkit"

	"""
	cnvkit2HRD.pl $segments ${meta.id} $ploidyv > ${meta.id}.cnvkit.scarHRD.txt
	"""   
}

process ASCAT2SCARHRD {
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/segments", mode: 'copy', overwrite: true, pattern: '*.txt'
 	cpus 1
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(id), val(type), file(logr), file(baf) //val(ploidy)

	output:
		tuple val(group), val(id), val(type), val("ascat"), file("${id}.ascat.scarHRD.txt"), emit: scarHRD_segments

	script:
		//ploidyv = ploidy.getText().trim()
		ploidyv = "NA"

	"""
	ascat2HRD.pl $logr $baf $id $ploidyv > ${id}.ascat.scarHRD.txt
	"""   
}

process ASCAT2OVAHRDSCAR {
	publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/segments", mode: 'copy', overwrite: true, pattern: '*.txt'
 	cpus 1
	time '1h'
	tag "$id"

	input:
		tuple val(group), val(id), val(type), file(logr), file(baf)

	output:
		tuple val(group), val(id), val(type), val("ascat"), file("${id}.ascat.ovaHRDscar.txt"), emit: ovaHRDscar_segments

	"""
	ascat2HRD.pl $logr $baf $id > ${id}.ascat.ovaHRDscar.txt
	"""   
}


process SCARHRD {
	publishDir "${params.outdir}/${params.subdir}/scarHRD/", mode: 'copy', overwrite: true, pattern: '*.txt'
	cpus 1
	time '1h'
	tag "${meta.id}"
	container = '/fs1/resources/containers/scarHRD.sif'

	input:
		tuple val(group), val(meta), val(sc), file(segments)

	output:
		tuple val(group), file("${meta.id}_${sc}_scarHRD_results.txt"), emit: scarHRD_score

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/Run_scarHRD.R $segments
	mv ${meta.id}_HRDresults.txt ${meta.id}_${sc}_scarHRD_results.txt
    """

}

process OVAHRDSCAR {
	publishDir "${params.outdir}/${params.subdir}/ovaHRDscar/", mode: 'copy', overwrite: true, pattern: '*.txt'
	cpus 1
	time '1h'
	tag "$id"
	container = '/fs1/resources/containers/ovaHRDscar.sif'

	input:
		tuple val(group), val(id), val(type), val(sc), file(segments)

	output:
		tuple val(group), val(id), val(type), file("${group}_${sc}_ovaHRDscar_results.txt"), emit: ovaHRDscar_score

	"""
	Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/Run_ovarHRDscar.R $segments > ${group}_${sc}_ovaHRDscar_results.txt
    """
}