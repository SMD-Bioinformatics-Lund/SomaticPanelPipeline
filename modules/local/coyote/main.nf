process COYOTE {
	publishDir "${params.crondir}/coyote", mode: 'copy', overwrite: true
	cpus 1
	time '10m'
	tag "$group"

	input:
		tuple val(group), val(meta), file(vcf), val(lowcov_type), file(lowcov), file(segments)

	output:
		tuple val(group), file("${process_group}.coyote"), emit: coyote_import

	when:
		!params.noupload

	script:
		process_group = group
		tumor_idx = 0
		tumor_idx_lowcov = 0
		if( meta.id.size() >= 2 ) {
			process_group = group + 'p'
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		}


		"""
		echo "import_myeloid_to_coyote_vep_gms.pl --group $params.coyote_group \\
			--vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
			--clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
			--lowcov /access/${params.subdir}/QC/${lowcov} \\
			--build 38 \\
			--gens ${group} \\
			--subpanel ${meta.diagnosis[tumor_idx]} \\
			--clarity-pool-id ${meta.clarity_pool_id[tumor_idx]}" > ${process_group}.coyote \\
			--cnv $segments \\
			--purity ${meta.purity[tumor_idx]}
		"""
	stub:
		process_group = group
		tumor_idx = 0
		tumor_idx_lowcov = 0
		if( meta.id.size() >= 2 ) {
			process_group = group + 'p'
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		}


		"""
		echo "import_myeloid_to_coyote_vep_gms.pl --group $params.coyote_group \\
			--vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
			--clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
			--lowcov /access/${params.subdir}/QC/${lowcov} \\
			--build 38 \\
			--gens ${group} \\
			--subpanel ${meta.diagnosis[tumor_idx]} \\
			--clarity-pool-id ${meta.clarity_pool_id[tumor_idx]}" > ${process_group}.coyote \\
			--cnv $segments \\
			--purity ${meta.purity[tumor_idx]}
		"""
}