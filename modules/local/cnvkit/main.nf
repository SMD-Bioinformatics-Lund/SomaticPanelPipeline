process CNVKIT {
	publishDir "${params.outdir}/${params.subdir}/cnvkit", mode: 'copy', overwrite: true, pattern: '*.cnvkit'
	publishDir "${params.outdir}/${params.subdir}/cnvkit/HRD", mode: 'copy', overwrite: true, pattern: '*.HRD'
	publishDir "${params.outdir}/${params.subdir}/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		tuple val(gr), val(id), val(type), file(bam), file(bai), file(bqsr), val(vc), file(vcf) 

	output:
		tuple val(gr), val(id), val(type), file("${gr}.${id}_logr_ballele.cnvkit"), emit: baflogr
		tuple file("${id}.HRD"), emit: cnvkit_hrd
		tuple file("${gr}.${id}.cnvkit_overview.png")

	when:
		params.cnvkit

	script:
		freebayes_idx = vc.findIndexOf{ it == 'freebayes' }

	"""
	set +eu
	source activate py2
	set -eu

	cnvkit.py batch $bam -r $params.cnvkit_reference -d results/
	cnvkit.py call results/*sort.cns -v $vcf -o ${gr}.${id}.call.cns	
	cnvkit.py scatter -s results/*sort.cn{s,r} -o ${gr}.${id}.cnvkit_overview.png -v ${vcf[freebayes_idx]} -i $id
	cp results/*sort.cnr ${gr}.${id}.cnr
	cp results/*sort.cns ${gr}.${id}.cns
	cnvkit.py export nexus-ogt -o 11-FF-HRD-GMSSTv1-0-210203.nexus -o ${gr}.${id}_logr_ballele.cnvkit ${gr}.${id}.cnr ${vcf[freebayes_idx]}
	generate_gens_data_from_cnvkit.pl ${gr}.${id}.cnr $vcf $id
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.cns tmp > ${id}.HRD
	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz" > ${id}.gens
	"""
}