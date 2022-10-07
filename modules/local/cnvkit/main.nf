process CNVKIT {
	publishDir "${params.outdir}/${params.subdir}/cnvkit", mode: 'copy', overwrite: true, pattern: '*.cnvkit'
	publishDir "${params.outdir}/${params.subdir}/cnvkit/segments", mode: 'copy', overwrite: true, pattern: '*.call.cns'
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
		call = "cnvkit.py call results/*sort.cns -v $vcf -o ${gr}.${id}.call.cns"
		if (purity) {
			call = call + "\ncnvkit.py call results/*sort.cns -v $vcf --purity $purity -o ${gr}.${id}.call.purity.cns"
		}

	"""
	set +eu
	source activate py2
	set -eu

	cnvkit.py batch $bam -r $params.cnvkit_reference -d results/
	$call
	cnvkit.py scatter -s results/*sort.cn{s,r} -o ${gr}.${id}.cnvkit_overview.png -v ${vcf} -i $id
	cp results/*sort.cnr ${gr}.${id}.cnr
	cp results/*sort.cns ${gr}.${id}.cns
	cnvkit.py export nexus-ogt -o 11-FF-HRD-GMSSTv1-0-210203.nexus -o ${gr}.${id}_logr_ballele.cnvkit ${gr}.${id}.cnr ${vcf}
	generate_gens_data_from_cnvkit.pl ${gr}.${id}.cnr $vcf $id
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.cns tmp > ${id}.HRD
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.purity.cns tmp > ${id}.purity.HRD
	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz" > ${id}.gens
	"""
}