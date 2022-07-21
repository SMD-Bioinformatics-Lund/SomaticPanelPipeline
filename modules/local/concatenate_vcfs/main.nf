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

process AGGREGATE_VCFS {
	cpus 1
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	time '20m'
	tag "$group"

	input:
		tuple val(group), val(vc), file(vcfs), val(id), val(type), val(FFPE)

	output:
		tuple val(group), file("${group}.agg.vcf"), emit: vcf_done

	script:
		sample_order = id[0]
		if( id.size() >= 2 ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			sample_order = id[tumor_idx]+","+id[normal_idx]
		}
		if (params.single_cnvcaller) {
			"""
			aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.vcf
			"""
		}
		else {
			"""
			aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.tmp.vcf
			vcf-concat ${group}.agg.tmp.vcf $cnvs | vcf-sort -c > ${group}.agg.vcf
			"""
		}

}