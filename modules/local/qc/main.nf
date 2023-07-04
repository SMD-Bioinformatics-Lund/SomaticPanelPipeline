process QC_TO_CDM {
	cpus 1
	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: 'true'
	tag "${meta.id}"
	time '10m'
	memory '50 MB'

	input:
		tuple val(group), val(meta), file(qc)

	output:
		tuple val(group), val(meta), file("${meta.id}.cdm"), emit: cdm_done

	when:
		!params.noupload

	script:
	"""
	echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${meta.id}.cdm
	"""
	stub:
	"""
	echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${meta.id}.cdm
	"""
}

process LOWCOV {
	cpus 1
	memory '5 GB'
	publishDir "${params.outdir}/${params.subdir}/QC", mode: 'copy', overwrite: 'true'
	time '1h'
	tag "${meta.id}"

	input:
		tuple val(group), val(meta), file(bam), file(bai), file(dedup) //from bam_lowcov

	output:
		tuple val(group), val(meta.type), file("${meta.id}.lowcov.bed"), emit: lowcov_regions

	script:
	"""
    source activate sambamba
	panel_depth.pl $bam ${params.regions_proteincoding} > lowcov.bed
	overlapping_genes.pl lowcov.bed ${params.gene_regions} > ${meta.id}.lowcov.bed
	"""
	stub:
	"""
	touch ${meta.id}.lowcov.bed
	"""
}

process QC_VALUES {
	time '2m'
	memory '50 MB'
	tag "${meta.id}"

	input:
		tuple val(group), val(meta), val(qc)

	output:
		tuple val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV), emit: qc_melt_val
	
	script:
		// Collect qc-data if possible from normal sample, if only tumor; tumor
        def ins_dev
        def coverage
        def ins_size
        qc.readLines().each{
			if (it =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
				ins_dev = it =~ /\"(ins_size_dev)\" : \"(\S+)\"/
			}
			if (it =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
				coverage = it =~ /\"(mean_coverage)\" : \"(\S+)\"/
			}
			if (it =~ /\"(ins_size)\" : \"(\S+)\"/) {
				ins_size = it =~ /\"(ins_size)\" : \"(\S+)\"/
			}
		}
		INS_SIZE = ins_size[0][2]
		MEAN_DEPTH = coverage[0][2]
		COV_DEV = ins_dev[0][2]
		"""
		echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
		"""
	stub:
		// Collect qc-data if possible from normal sample, if only tumor; tumor
        def ins_dev
        def coverage
        def ins_size
		INS_SIZE = 1
		MEAN_DEPTH = 1
		COV_DEV = 1
		"""
		echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
		"""
}
