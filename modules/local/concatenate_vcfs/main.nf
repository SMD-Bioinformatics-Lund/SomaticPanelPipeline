process CONCATENATE_VCFS {
	cpus 1
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	time '20m'    
	tag "$group"

	input:
		tuple val(vc), val(group), file(vcfs)

	output:
		tuple val(group), val(vc), file("${group}_${vc}.vcf.gz"), emit: concatenated_vcfs

	"""
	vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt normalize ${vc}.decomposed.vcf.gz -r $params.genome_file | vt uniq - -o ${group}_${vc}.vcf.gz
	"""
}