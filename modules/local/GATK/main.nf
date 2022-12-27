process GATKCOV_BAF {
	cpus 10
	memory '64 GB'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

	input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

	output:
        tuple val(group), val(meta), file("${meta.id}.allelicCounts.tsv"), emit: gatk_baf
    
	script:
		
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx50g" CollectAllelicCounts \\
			-L $params.GATK_GNOMAD \\
			-I $bam \\
			-R $params.genome_file \\
			-O ${meta.id}.allelicCounts.tsv
		"""
	stub:
		"""
		touch ${meta.id}.allelicCounts.tsv
		"""


}

process GATKCOV_COUNT {
	cpus 10
	memory '64 GB'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true 

	input:
		tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

	output:
		tuple val(group), val(meta), file("${meta.id}.standardizedCR.tsv"), file("${meta.id}.denoisedCR.tsv"), emit: gatk_count

	script:

		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk CollectReadCounts \\
			-I ${bam} -L $params.gatk_intervals \\
			--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5
		gatk --java-options '-Xmx50g' DenoiseReadCounts \\
			-I ${bam}.hdf5 --count-panel-of-normals $params.GATK_pon \\
			--standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${meta.id}.denoisedCR.tsv
		gatk PlotDenoisedCopyRatios \\
			--standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${meta.id}.denoisedCR.tsv \\
			--sequence-dictionary $params.GENOMEDICT \\
			--minimum-contig-length 46709983 --output . --output-prefix ${meta.id}
		"""
	stub:
		"""
		touch ${meta.id}.standardizedCR.tsv ${meta.id}.denoisedCR.tsv
		"""	
}

process GATKCOV_CALL {
	publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true 
	cpus 10
	memory '64 GB'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

	input:
		tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

	output:
		tuple val(group), val(meta), file("${meta.id}.called.seg"), emit: gatcov_called
		tuple val(group), val(meta), file("${meta.id}.modeled.png"), emit: gatcov_plot

	script:
		// tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T'  }
		// normal_idx = type.findIndexOf{ it == 'normal' || it == 'N'  }

	// --normal-allelic-counts ${id[normal_idx]}.allelicCounts.tsv \\
	// add to modelsegments when paired? //
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx40g" ModelSegments \\
			--denoised-copy-ratios $denoised \\
			--allelic-counts $allele \\
			--minimum-total-allele-count-normal 20 \\
			--output . \\
			--output-prefix ${meta.id}
		gatk CallCopyRatioSegments \\
			--input ${meta.id}.cr.seg \\
			--output ${meta.id}.called.seg
		gatk PlotModeledSegments \\
			--denoised-copy-ratios $denoised \\
			--allelic-counts ${meta.id}.hets.tsv \\
			--segments ${meta.id}.modelFinal.seg \\
			--sequence-dictionary $params.GENOMEDICT \\
			--minimum-contig-length 46709983 \\
			--output . \\
			--output-prefix ${meta.id}
		"""
	stub:
		"""
		touch ${meta.id}.called.seg
		touch ${meta.id}.modeled.png
		"""
}

