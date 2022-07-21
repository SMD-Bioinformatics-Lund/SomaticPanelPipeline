process MSISENSOR {
        cpus 2
        memory '1GB'
        publishDir "${params.outdir}/${params.subdir}/msi", mode: 'copy', overwrite: true
        time '20m'
        tag "$group"
		container = "/fs1/resources/containers/msisensor-pro-1.2.0.sif"

		when:
			params.msi

        input:
			tuple val(group), val(id), val(type), file(bams), file(bais), file(bqsr)
			
        output:
            tuple val(group), file("${group}.msi*"), emit: msi_score

		script:

			if( id.size() >= 2 ) {

				tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
				normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

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
			else if( id.size() == 1 ) {
				"""
				msisensor-pro pro \\
					-d /fs1/resources/validation-samples/solid/msisensorbaseline/msisensor_reference_hg38.list_baseline \\
					-t $bams \\
					-o ${group}.msi_single -c 50 -b ${task.cpus}
				"""
			}

}