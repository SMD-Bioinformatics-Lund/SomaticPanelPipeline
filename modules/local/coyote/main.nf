process COYOTE {
	publishDir "${params.crondir}/coyote", mode: 'copy', overwrite: true
	cpus 1
	time '10m'
	tag "$group"

	input:
		tuple val(group), file(vcf), val(id), val(type), val(lims_id), val(pool_id), val(diagnosis), val(id_lowcov), val(type_lowcov), file(lowcov)

	output:
		tuple val(group), file("${process_group}.coyote"), emit: coyote_import

	when:
		!params.noupload

	script:
		process_group = group
		tumor_idx = 0
		tumor_idx_lowcov = 0
		if( id.size() >= 2 ) {
			process_group = group + 'p'
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			tumor_idx_lowcov = type_lowcov.findIndexOf{ it == 'tumor' || it == 'T' }
		}


	"""
	echo "import_myeloid_to_coyote_vep_gms.pl --group $params.coyote_group \\
		--vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
		--clarity-sample-id ${lims_id[tumor_idx]} \\
		--lowcov /access/${params.subdir}/QC/${lowcov[tumor_idx_lowcov]} \\
		--build 38 \\
		--gens ${group} \\
		--subpanel ${diagnosis[tumor_idx]} \\
		--clarity-pool-id ${pool_id[tumor_idx]}" > ${process_group}.coyote
	"""
}