process VARDICT {
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
		each file(bed)

	output:
		tuple val("vardict"), val(group), file("vardict_${bed}.vcf"), emit: vcfparts_vardict

	when:
		params.vardict
    
	script:
		if( meta.id.size() >= 2 ) {

			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			vardict-java -G $params.genome_file -f 0.01 -N ${meta.id[tumor_idx]} -b "${bams[tumor_idx]}|${bams[normal_idx]}" -c 1 -S 2 -E 3 -g 4 -U $bed \\
			| testsomatic.R | var2vcf_paired.pl -N "${meta.id[tumor_idx]}|${meta.id[normal_idx]}" -f 0.01 > vardict_${bed}.vcf.raw

			filter_vardict_somatic.pl vardict_${bed}.vcf.raw ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > vardict_${bed}.vcf
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			vardict-java -G $params.genome_file -f 0.03 -N ${meta.id[0]} -b ${bams[0]} -c 1 -S 2 -E 3 -g 4 -U $bed | teststrandbias.R | var2vcf_valid.pl -N ${meta.id[0]} -E -f 0.01 > vardict_${bed}.vcf.raw
			filter_vardict_unpaired.pl vardict_${bed}.vcf.raw > vardict_${bed}.vcf
			"""
		}
	stub:
	if( meta.id.size() >= 2 ) {
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
		"""
		echo tumor:${bams[tumor_idx]} ${meta.id[tumor_idx]} normal:${bams[normal_idx]} ${meta.id[normal_idx]}
		touch vardict_${bed}.vcf
		"""
	}
	else {
		"""
		echo tumor:$bams
		touch vardict_${bed}.vcf
		"""
	}

}