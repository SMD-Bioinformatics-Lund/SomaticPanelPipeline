process CNVKIT {
	publishDir "${params.outdir}/${params.subdir}/cnvkit", mode: 'copy', overwrite: true, pattern: '*.cnvkit'
	publishDir "${params.outdir}/${params.subdir}/cnvkit/segments", mode: 'copy', overwrite: true, pattern: '*.call.cns'
	publishDir "${params.outdir}/${params.subdir}/cnvkit/HRD", mode: 'copy', overwrite: true, pattern: '*.HRD'
	publishDir "${params.outdir}/${params.subdir}/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	//scratch true
	//stageInMode 'copy'
	//stageOutMode 'copy'
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		//tuple val(gr), val(id), val(type), file(bam), file(bai), file(bqsr), val(purity), val(vc), file(vcf)
		tuple val(gr), val(id), val(type), file(bam), file(bai), file(bqsr), val(purity), file(vcf), file(tbi)

	output:
		tuple val(gr), val(id), val(type), file("${gr}.${id}_logr_ballele.cnvkit"), emit: baflogr
		tuple val(gr), val(id), val(type), file("${id}.HRD"), file("${id}.purity.HRD"), emit: cnvkit_hrd
		tuple val(gr), val(id), val(type), file("${gr}.${id}.cnvkit_overview.png")
		tuple val(gr), val(id), val(type), file("${gr}.${id}.call.cns"), emit: cnvkitsegment
		tuple val(gr), val(id), val(type), file("${gr}.${id}.call.purity.cns"), emit: cnvkitsegment_purity

	when:
		params.cnvkit

	script:
		//freebayes_idx = vc.findIndexOf{ it == 'freebayes' }
		call = "cnvkit.py call results_backbone/*sort.cns -v $vcf -o ${gr}.${id}.call.cns"
		if (purity) {
			call = call + "\ncnvkit.py call results_backbone/*sort.cns -v $vcf --purity $purity -o ${gr}.${id}.call.purity.cns"
		}

	"""
	set +eu
	source activate py2
	set -eu

	cnvkit.py batch $bam -r $params.cnvkit_reference -d results/
	cnvkit.py batch $bam -r /fs2/viktor/CNVkit/ref_SOLID_v3/solidv3_17normals_exons.cnn -d results_exons/
	cnvkit.py batch $bam -r /fs2/viktor/CNVkit/ref_SOLID_v3/solidv3_17normals_backbone.cnn -d results_backbone/
	$call
	cnvkit.py scatter -s results_backbone/*sort.cn{s,r} -o ${gr}.${id}.cnvkit_overview.png -v ${vcf} -i $id
	cp results_backbone/*sort.cnr ${gr}.${id}.cnr
	cp results_backbone/*sort.cns ${gr}.${id}.cns
	cnvkit.py export nexus-ogt -o 11-FF-HRD-GMSSTv1-0-210203.nexus -o ${gr}.${id}_logr_ballele.cnvkit ${gr}.${id}.cnr ${vcf}
	generate_gens_data_from_cnvkit.pl ${gr}.${id}.cnr $vcf $id
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.cns tmp > ${id}.HRD
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.purity.cns tmp > ${id}.purity.HRD
	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz" > ${id}.gens
	"""
}
process CNVKIT_BATCH {
	publishDir "${params.outdir}/${params.subdir}/cnvkit", mode: 'copy', overwrite: true
	cpus 2
	time '1h'
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		tuple val(group), val(meta), file(bam), file(bai), file(bqsr)
		val(reference)

	output:
		tuple val(group), val(meta), file("${group}.${meta.id}.cns"), emit: cnvkit_cns
		tuple val(group), val(meta), file("${group}.${meta.id}.cnr"), emit: cnvkit_cnr
		

	when:
		params.cnvkit

	script:
		"""
		cnvkit.py batch $bam -r $reference -d results/
		cp results/*sort.cnr ${group}.${meta.id}.cnr
		cp results/*sort.cns ${group}.${meta.id}.cns
		"""
	stub:
		"""
		touch ${group}.${meta.id}.cns ${group}.${meta.id}.cnr
		echo $reference
		"""
}

process CNVKIT_PLOT {
	publishDir "${params.outdir}/${params.subdir}/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		tuple val(group), val(meta), file(cns), file(cnr), file(vcf), file(tbi)

	output:
		tuple val(group), val(meta), file("${group}.${meta.id}.cnvkit_overview.png")

	when:
		params.cnvkit

	script:
		"""
		echo $cnr $cns
		cnvkit.py scatter -s *.cn{s,r} -o ${group}.${meta.id}.cnvkit_overview.png -v ${vcf} -i ${meta.id}
		"""
	stub:
		"""
		touch ${group}.${meta.id}.cnvkit_overview.png
		"""
}

process CNVKIT_GENS {
	publishDir "${params.outdir}/${params.subdir}/gens", mode: 'copy', overwrite: true, pattern: '*.bed.gz'
	cpus 1
	time '1h'
	tag "${meta.id}"
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		tuple val(group), val(meta), file(cnr), file(vcf), file(tbi)

	output:
		tuple val(group), file("${meta.id}.baf.bed.gz"), file("${meta.id}.cov.bed.gz"), emit: cnvkit_gens
		tuple val(group), file("${meta.id}.gens"), emit: cnvkit_gens_loaddb
		
	when:
		params.cnvkit

	script:
		"""
		set +eu
		source activate py2
		set -eu
		generate_gens_data_from_cnvkit.pl ${group}.${meta.id}.cnr $vcf ${meta.id}
		echo "gens load sample --sample-id ${meta.id} --genome-build 38 --baf ${params.gens_accessdir}/${meta.id}.baf.bed.gz --coverage ${params.gens_accessdir}/${meta.id}.cov.bed.gz" > ${meta.id}.gens
		"""
	stub:
		"""
		touch ${meta.id}.baf.bed.gz ${meta.id}.cov.bed.gz ${meta.id}.gens
		"""
}

process CNVKIT_CALL {
	publishDir "${params.outdir}/${params.subdir}/cnvkit/segments/", mode: 'copy', overwrite: true, pattern: '*call*.cns'
	cpus 1
	time '1h'
	tag "${meta.id}"
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		tuple val(group), val(meta), file(cns), file(cnr), file(vcf), file(tbi)

	output:
		tuple val(group), val(meta), file("${group}.${meta.id}.call*.cns"), emit: cnvkitsegment
		tuple val(group), val(meta), file("${group}.${meta.id}_logr_ballele.cnvkit"), emit: cnvkit_baflogr
		
	when:
		params.cnvkit

	script:
		call = "cnvkit.py call $cns -v $vcf -o ${group}.${meta.id}.call.cns"
		if (meta.purity) {
			call = "cnvkit.py call $cns -v $vcf --purity ${meta.purity} -o ${group}.${meta.id}.call.purity.cns"
		}
		"""
		set +eu
		source activate py2
		set -eu
		$call
		cnvkit.py export nexus-ogt -o ${group}.${meta.id}_logr_ballele.cnvkit ${cnr} ${vcf}
		"""
	stub:
		call = "cnvkit.py call $cns -v $vcf -o ${group}.${meta.id}.call.cns"
		if (meta.purity) {
			call = "cnvkit.py call $cns -v $vcf --purity ${meta.purity} -o ${group}.${meta.id}.call.purity.cns"
		}
		"""
		echo $call
		touch ${group}.${meta.id}.call.purity.cns ${group}.${meta.id}_logr_ballele.cnvkit
		"""
}