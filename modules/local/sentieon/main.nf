process BWA_UMI {
	publishDir "${params.outdir}/${params.subdir}/bam", mode: 'copy', overwrite: true, pattern: "*.bam*"
	cpus params.cpu_all
	memory '128 GB'
	time '2h'
	errorStrategy 'retry'
	maxErrors 5
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(meta), file(r1), file(r2)

	output:
		tuple val(group), val(meta), file("${out_umi}"), file("${out_umi}.bai"), emit: bam_umi
		tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), emit: bam_umi_markdup

	when:
		params.umi

	script:
		out_bam = meta.id+"."+meta.type+".bwa.sort.bam"
		out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"
		if (meta.sub) {
			submbp = params.sample_val / 1000000
			submbp = submbp + "M"
			out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
			out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
		}
		"""

		export skip_coord_end=true
		
		sentieon umi extract -d 3M2S+T,3M2S+T $r1 $r2 \\
		|sentieon bwa mem \\
			-R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:illumina" \\
			-t ${task.cpus} \\
			-p -C $params.genome_file - \\
		|tee -a noumi.sam \\
		|sentieon umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

		sentieon bwa mem \\
			-R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:illumina" \\
			-t ${task.cpus} \\
			-p -C $params.genome_file consensus.fastq.gz \\
		|sentieon util sort -i - \\
			-o ${meta.id}.${meta.type}.bwa.umi.sort.bam \\
			--sam2bam --umi_post_process

		sentieon util sort -i noumi.sam -o ${out_bam} --sam2bam
		rm noumi.sam

		touch dedup_metrics.txt
		"""
	stub:
		out_bam = meta.id+"."+meta.type+".bwa.sort.bam"
		out_umi = meta.id+"."+meta.type+".bwa.umi.sort.bam"
		if (meta.sub) {
			submbp = params.sample_val / 1000000
			submbp = submbp + "M"
			out_bam = meta.id+"."+meta.type+"."+submbp+".bwa.sort.bam"
			out_umi = meta.id+"."+meta.type+"."+submbp+".bwa.umi.sort.bam"
		}
		"""
		touch ${out_bam} ${out_bam}.bai
		touch ${out_umi} ${out_umi}.bai
		"""
}

process MARKDUP {
	publishDir "${params.outdir}/${params.subdir}/bam", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '64 GB'
	time '1h'
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"
	
	input:
		tuple val(group), val(meta), file(bam), file(bai)

	output:
		tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), emit: bam_bqsr
		tuple val(group), val(meta), file("${out_bam}"), file("${out_bam}.bai"), file("dedup_metrics.txt"), emit: bam_qc

	script:
		out_bam = meta.id+"."+meta.type+".dedup.bam"
		if (meta.sub) {
			submbp = params.sample_val / 1000000
			submbp = submbp + "M"
			out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
		}
		"""
		sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info score.gz
		sentieon driver -t ${task.cpus} -i $bam --algo Dedup --score_info score.gz --metrics dedup_metrics.txt $out_bam
		"""
	stub:
		out_bam = meta.id+"."+meta.type+".dedup.bam"
		if (meta.sub) {
			submbp = params.sample_val / 1000000
			submbp = submbp + "M"
			out_bam = meta.id+"."+meta.type+"."+submbp+".dedup.bam"
		}
		"""
		touch ${out_bam} ${out_bam}.bai dedup_metrics.txt
		"""
}

process BQSR_UMI {
	cpus params.cpu_some
	memory '16 GB'
	time '1h'
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(meta), file(bam), file(bai)

	output:
		tuple val(group), val(meta), file(bam), file(bai), file("${meta.id}.bqsr.table"), emit: bam_varcall
	when:
		params.umi
	script:
		"""
		sentieon driver -t ${task.cpus} -r $params.genome_file -i $bam --algo QualCal ${meta.id}.bqsr.table
		"""
	stub:
		"""
		touch ${meta.id}.bqsr.table
		"""
}

process SENTIEON_QC {
	cpus params.cpu_many
	memory '32 GB'
	publishDir "${params.outdir}/${params.subdir}/QC", mode: 'copy', overwrite: 'true', pattern: '*.QC*'
	time '1h'
	tag "${meta.id}"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(meta), file(bam), file(bai), file(dedup)

	output:
		tuple val(group), val(meta), file(bam), file(bai), file("${meta.id}_is_metrics.txt"), emit: dedup_bam
		tuple val(group), val(meta), file("${meta.id}_${meta.type}.QC"), emit: qc_cdm
		tuple val(group), val(meta), file("${meta.id}_${meta.type}.QC"), emit: qc_melt
		file("*.txt")

	script:
		"""
		sentieon driver \\
			--interval $params.regions_bed_qc -r $params.genome_file -t ${task.cpus} -i ${bam} \\
			--algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \\
			--algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt \\
			--algo InsertSizeMetricAlgo is_metrics.txt \\
			--algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
		sentieon driver \\
			-r $params.genome_file -t ${task.cpus} -i ${bam} \\
			--algo HsMetricAlgo --targets_list $params.interval_list_qc --baits_list $params.interval_list_qc hs_metrics.txt

		cp is_metrics.txt ${meta.id}_is_metrics.txt

		qc_sentieon.pl ${meta.id}_${meta.type} panel > ${meta.id}_${meta.type}.QC
		"""
	stub:
		"""
		touch ${meta.id}_is_metrics.txt
		touch ${meta.id}_${meta.type}.QC
		"""
}

