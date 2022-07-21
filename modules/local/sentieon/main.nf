process BWA_UMI {
	publishDir "${params.outdir}/${params.subdir}/bam", mode: 'copy', overwrite: true, pattern: "*.bam*"
	cpus params.cpu_all
	memory '128 GB'
	time '2h'
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(id), val(type), file(r1), file(r2)

	output:
		tuple val(group), val(id), val(type), file("${id}.${type}.bwa.umi.sort.bam"), file("${id}.${type}.bwa.umi.sort.bam.bai"), emit: bam_umi
		tuple val(group), val(id), val(type), file("${id}.${type}.bwa.sort.bam"), file("${id}.${type}.bwa.sort.bam.bai"), emit: bam_umi_markdup

	when:
		params.umi

	"""

	export skip_coord_end=true
	
	sentieon umi extract -d 3M2S+T,3M2S+T $r1 $r2 \\
	|sentieon bwa mem \\
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \\
		-t ${task.cpus} \\
		-p -C $params.genome_file - \\
	|tee -a noumi.sam \\
	|sentieon umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

	sentieon bwa mem \\
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \\
		-t ${task.cpus} \\
		-p -C $params.genome_file consensus.fastq.gz \\
	|sentieon util sort -i - \\
		-o ${id}.${type}.bwa.umi.sort.bam \\
		--sam2bam --umi_post_process

	sentieon util sort -i noumi.sam -o ${id}.${type}.bwa.sort.bam --sam2bam
	rm noumi.sam

	touch dedup_metrics.txt
	"""
}

process MARKDUP {
	publishDir "${params.outdir}/${params.subdir}/bam", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '64 GB'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    container = "/fs1/resources/containers/sentieon_202112.sif"
    
	input:
		tuple val(group), val(id), val(type), file(bam), file(bai)

	output:
		tuple val(group), val(id), val(type), file("${id}.${type}.dedup.bam"), file("${id}.${type}.dedup.bam.bai"), emit: bam_bqsr
		tuple val(group), val(id), val(type), file("${id}.${type}.dedup.bam"), file("${id}.${type}.dedup.bam.bai"), file("dedup_metrics.txt"), emit: bam_qc

	"""
	sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info score.gz
	sentieon driver -t ${task.cpus} -i $bam --algo Dedup --score_info score.gz --metrics dedup_metrics.txt ${id}.${type}.dedup.bam
	"""
}

process BQSR_UMI {
	cpus params.cpu_some
	memory '16 GB'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(id), val(type), file(bam), file(bai)

	output:
		tuple val(group), val(id), val(type), file(bam), file(bai), file("${id}.bqsr.table"), emit: bam_varcall
	when:
		params.umi

	"""
	sentieon driver -t ${task.cpus} -r $params.genome_file -i $bam --algo QualCal ${id}.bqsr.table
	"""
}

process SENTIEON_QC {
	cpus params.cpu_many
	memory '32 GB'
	publishDir "${params.outdir}/${params.subdir}/QC", mode: 'copy', overwrite: 'true', pattern: '*.QC*'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		tuple val(group), val(id), val(type), file(bam), file(bai), file(dedup)

	output:
		tuple val(group), val(id), val(type), file(bam), file(bai), file("${id}_is_metrics.txt"), emit: dedup_bam
		tuple val(id), val(type), file("${id}_${type}.QC"), emit: qc_cdm
		tuple val(group), val(id), val(type), file("${id}_${type}.QC"), emit: qc_melt
		file("*.txt")

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

	cp is_metrics.txt ${id}_is_metrics.txt

	qc_sentieon.pl ${id}_${type} panel > ${id}_${type}.QC
	"""
}