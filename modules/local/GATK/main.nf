process GATKCOV_BAF {
	cpus 10
	memory '64 GB'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

	input:
        tuple val(group), val(id), val(type), file(bam), file(bai), file(bqsr)

	output:
        tuple val(group), val(id), val(type), file("${id}.allelicCounts.tsv"), emit: gatk_baf
    
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
		-O ${id}.allelicCounts.tsv
	"""

}

process GATKCOV_COUNT {
	cpus 10
	memory '64 GB'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true 

	input:
		tuple val(group), val(id), val(type), file(bam), file(bai), file(bqsr)

	output:
		tuple val(group), val(id), val(type), file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv"), emit: gatk_count

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
		--standardized-copy-ratios ${id}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id}.denoisedCR.tsv
	gatk PlotDenoisedCopyRatios \\
		--standardized-copy-ratios ${id}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id}.denoisedCR.tsv \\
		--sequence-dictionary $params.GENOMEDICT \\
		--minimum-contig-length 46709983 --output . --output-prefix ${id}
    """
}

process GATKCOV_CALL {
	publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true 
	cpus 10
	memory '64 GB'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

	input:
		tuple val(group), val(id), val(type), file(allele), file(stdCR), file(denoised)

	output:
		

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
		--output-prefix ${id}
	gatk CallCopyRatioSegments \\
		--input ${id}.cr.seg \\
		--output ${id}.called.seg
	gatk PlotModeledSegments \\
		--denoised-copy-ratios $denoised \\
		--allelic-counts ${id}.hets.tsv \\
		--segments ${id}.modelFinal.seg \\
		--sequence-dictionary $params.GENOMEDICT \\
		--minimum-contig-length 46709983 \\
		--output . \\
		--output-prefix ${id}
	"""
}

