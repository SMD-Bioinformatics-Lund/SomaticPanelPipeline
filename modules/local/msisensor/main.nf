process MSISENSOR {
	cpus 2
	memory '1GB'
	publishDir "${params.outdir}/${params.subdir}/msi", mode: 'copy', overwrite: true
	time '20m'
	tag "$group"
	container = "/fs1/resources/containers/msisensor-pro-1.2.0.sif"
	errorStrategy 'ignore'

	when:
		params.msi

	input:
		tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
		
	output:
		tuple val(group), file("${group}.msi_single"), emit: msi_score
		tuple val(group), file("${group}.msi_paired"), optional:true, emit: msi_score_paired

	script:

		if( meta.id.size() >= 2 ) {

			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			msisensor-pro msi \\
				-d /fs1/resources/validation-samples/solid/msisensorbaseline/msisensor_reference_hg38.list \\
				-n ${bams[normal_idx]} \\
				-t ${bams[tumor_idx]} \\
				-o ${group}.msi_paired -c 50 -b ${task.cpus}
			
			msisensor-pro pro \\
				-d /fs1/resources/validation-samples/solid/msisensorbaseline/msisensor_reference_hg38.list_baseline \\
				-t ${bams[tumor_idx]} \\
				-o ${group}.msi_single -c 50 -b ${task.cpus}
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			msisensor-pro pro \\
				-d /fs1/resources/validation-samples/solid/msisensorbaseline/msisensor_reference_hg38.list_baseline \\
				-t $bams \\
				-o ${group}.msi_single -c 50 -b ${task.cpus}
			"""
		}
	stub:
		if( meta.id.size() >= 2 ) {

			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			touch ${group}.msi_single
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			touch ${group}.msi_single
			"""
		}
}