process CONCATENATE_VCFS {
	label "process_single"
	tag "$group"

	input:
		tuple val(vc), val(group), file(vcfs)

	output:
		tuple val(group), val(vc), file("${group}_${vc}.vcf.gz"), emit: concatenated_vcfs

	script:
	"""
	vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt normalize ${vc}.decomposed.vcf.gz -r $params.genome_file | vt uniq - -o ${group}_${vc}.vcf.gz
	"""
	stub:
	"""
	touch ${group}_${vc}.vcf.gz
	"""
}

process AGGREGATE_VCFS {
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(vc), file(vcfs), val(meta)

	output:
		tuple val(group), val(meta), file("${group}.agg.vcf"), emit: vcf_concat

	script:
		sample_order = meta.id[0]
		if( meta.id.size() >= 2 ) {
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
		}

		"""
		aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.vcf
		"""


	stub:
	if( meta.id.size() >= 2 ) {
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
		sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
		"""
		echo tumor:${meta.id[tumor_idx]} normal:${meta.id[normal_idx]}
		touch ${group}.agg.vcf
		"""
	}
	else {
		"""
		touch ${group}.agg.vcf
		"""
	}

}
