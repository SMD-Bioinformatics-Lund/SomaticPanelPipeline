#!/usr/bin/env nextflow

genome_file = file(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

csv = file(params.csv)
println(csv)



workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}

// Print commit-version of active deployment
file(params.git)
    .readLines()
    .each { println "git commit-hash: "+it }
// Print active container
container = file(params.container).toRealPath()
println("container: "+container)

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, file(row.read1), file(row.read2)) }
    .into { fastq_umi; fastq_noumi; meta_nocnv }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, (row.containsKey("ffpe") ? row.ffpe : false)) }
    .into { meta_aggregate; meta_germline; meta_pon; meta_cnvkit; meta_melt; meta_cnvplot }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.type, row.clarity_sample_id, row.clarity_pool_id) }
    .set { meta_coyote }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, row.read1, row.read2) }
    .set{ meta_qc }

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, row.clarity_sample_id, row.clarity_pool_id) }
    .set { meta_const }



// Split bed file in to smaller parts to be used for parallel variant calling
Channel
    .fromPath("${params.regions_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.regions_bed}" }
    .splitText( by: 200, file: 'bedpart.bed' )
    .into { beds_mutect; beds_freebayes; beds_tnscope; beds_vardict }

Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.into{ gatk_ref; gatk_postprocess }

process bwa_umi {
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: true, pattern: "*.bam*"
	cpus params.cpu_all
	memory '128 GB'
	time '2h'
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, type, file(r1), file(r2) from fastq_umi

	output:
		set group, id, type, file("${id}.${type}.bwa.umi.sort.bam"), file("${id}.${type}.bwa.umi.sort.bam.bai") into bam_umi_bqsr, bam_umi_confirm
		set group, id, type, file("${id}.${type}.bwa.sort.bam"), file("${id}.${type}.bwa.sort.bam.bai") into bam_umi_markdup

	when:
		params.umi

	"""

	export skip_coord_end=true
	
	sentieon umi extract -d 3M2S+T,3M2S+T $r1 $r2 \\
	|sentieon bwa mem \\
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \\
		-t ${task.cpus} \\
		-p -C $genome_file - \\
	|tee -a noumi.sam \\
	|sentieon umi consensus --copy_tags XR,RX,MI,XZ -o consensus.fastq.gz

	sentieon bwa mem \\
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \\
		-t ${task.cpus} \\
		-p -C $genome_file consensus.fastq.gz \\
	|sentieon util sort -i - \\
		-o ${id}.${type}.bwa.umi.sort.bam \\
		--sam2bam --umi_post_process

	sentieon util sort -i noumi.sam -o ${id}.${type}.bwa.sort.bam --sam2bam
	rm noumi.sam

	touch dedup_metrics.txt
	"""
}


process bwa_align {
	cpus params.cpu_all
	memory '64 GB'
	time '2h'
	tag "$id"
	    
	input: 
		set group, id, type, file(r1), file(r2) from fastq_noumi

	output:
		set group, id, type, file("${id}.${type}.bwa.sort.bam"), file("${id}.${type}.bwa.sort.bam.bai") into bam_markdup

	when:
		!params.umi

	script:

		if( params.sentieon_bwa ) {
			"""
			sentieon bwa mem -M -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -t ${task.cpus} $genome_file $r1 $r2 \\
			| sentieon util sort -r $genome_file -o ${id}.${type}.bwa.sort.bam -t ${task.cpus} --sam2bam -i -
			"""
		}

		else {
			"""
			bwa mem -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -M -t ${task.cpus} $genome_file $r1 $r2 \\
			| samtools view -Sb - \\
			| samtools sort -o ${id}.${type}.bwa.sort.bam -

			samtools index ${id}.${type}.bwa.sort.bam
			"""
		}
}


process markdup {
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '64 GB'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    
	input:
		set group, id, type, file(bam), file(bai) from bam_markdup.mix(bam_umi_markdup)

	output:
		set group, id, type, file("${id}.${type}.dedup.bam"), file("${id}.${type}.dedup.bam.bai") into bam_bqsr
		set group, id, type, file("${id}.${type}.dedup.bam"), file("${id}.${type}.dedup.bam.bai"), file("dedup_metrics.txt") into bam_qc, bam_bqsr2, bam_lowcov

	"""
	sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info score.gz
	sentieon driver -t ${task.cpus} -i $bam --algo Dedup --score_info score.gz --metrics dedup_metrics.txt ${id}.${type}.dedup.bam
	"""
}

// FIXME: Temporarily broke the non-UMI track since bam_umi_bqsr
//        and bam_bqsr collide here for UMI track. Figure out how
//        to use only bam_umi_bqsr when params.umi==true
process bqsr_umi {
	cpus params.cpu_some
	memory '16 GB'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, type, file(bam), file(bai) from bam_umi_bqsr

	output:
		set group, id, type, file(bam), file(bai), file("${id}.bqsr.table") into bam_freebayes, bam_vardict, bam_tnscope, bam_cnvkit, bam_varli, bam_gatk
	when:
		params.umi

	"""
	sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo QualCal ${id}.bqsr.table
	"""
}

process bqsr_to_constitutional {
	cpus params.cpu_some
	memory '16 GB'
	time '1h'
	tag "$id"
	publishDir "${OUTDIR}/bqsr", mode: 'copy', overwrite: true, pattern: '*.bqsr*'
	publishDir "${OUTDIR}/csv", mode: 'copy', overwrite: true, pattern: '*.csv*'
	//scratch true
	//stageInMode 'copy'
	//stageOutMode 'copy'

	when:
		mode == "neverever"

	input:
		set group, id, type, file(bam), file(bai), file(dedup), cid, poolid from bam_bqsr2.join(meta_const, by: [0,1,2]).filter{ item -> item[2] == 'N'}

	output:
		set group, id, type, file(bam), file(bai), file("${id}.const.bqsr") into input_const
		file("${id}.const.csv") into csv_const

	"""
	sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo QualCal ${id}.const.bqsr
	cat $params.const_csv_template > ${id}.const.csv
	echo $cid,$id,proband,oncov1-0-test,M,ovarian-normal,affected,$id,,,$poolid,illumina,${OUTDIR}/bam/$bam,${OUTDIR}/bqsr/${id}.const.bqsr,,screening >> ${id}.const.csv
	/fs1/bjorn/bnf-scripts/start_nextflow_analysis.pl ${id}.const.csv
	"""
	
}

process sentieon_qc {
	cpus params.cpu_many
	memory '32 GB'
	publishDir "${OUTDIR}/QC", mode: 'copy', overwrite: 'true', pattern: '*.QC*'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, type, file(bam), file(bai), file(dedup) from bam_qc

	output:
		set group, id, type, file(bam), file(bai), file("${id}_is_metrics.txt") into all_pindel, bam_manta, bam_melt, bam_delly
		set id, type, file("${id}_${type}.QC") into qc_cdm
		set group, id, type, file("${id}_${type}.QC") into qc_melt
		file("*.txt")

	"""
	sentieon driver \\
		--interval $params.regions_bed_qc -r $genome_file -t ${task.cpus} -i ${bam} \\
		--algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \\
		--algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt \\
		--algo InsertSizeMetricAlgo is_metrics.txt \\
		--algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
	sentieon driver \\
		-r $genome_file -t ${task.cpus} -i ${bam} \\
		--algo HsMetricAlgo --targets_list $params.interval_list_qc --baits_list $params.interval_list_qc hs_metrics.txt

	cp is_metrics.txt ${id}_is_metrics.txt

	qc_sentieon.pl ${id}_${type} panel > ${id}_${type}.QC
	"""
}

process lowcov {
	cpus 1
	memory '5 GB'
	publishDir "${OUTDIR}/QC", mode: 'copy', overwrite: 'true'
	time '1h'
	tag "$id"

	input:
		set group, id, type, file(bam), file(bai), file(dedup) from bam_lowcov


	output:
		set group, type, file("${id}.lowcov.bed") into lowcov_coyote

	"""
    source activate sambamba
	panel_depth.pl $bam $params.regions_proteincoding > lowcov.bed
	overlapping_genes.pl lowcov.bed $params.gene_regions > ${id}.lowcov.bed
	"""
}

// Load QC data into CDM (via middleman)
process qc_to_cdm {
	cpus 1
	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	time '10m'
	memory '50 MB'

	input:
		set id, type, file(qc), r1, r2 from qc_cdm.join(meta_qc)

	output:
		file("${id}.cdm") into cdm_done

	when:
		!params.noupload

	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")

	"""
	echo "--run-folder $rundir --sample-id $id --assay $params.cdm --qc ${OUTDIR}/QC/${id}_${type}.QC" > ${id}.cdm
	"""
}

process qc_values {
	tag "$id"
	time '2m'
	memory '50 MB'
	tag "$id"

	input:
		set group, id, type, qc from qc_melt

	output:
		set group, id, type, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_melt_val
		set group, id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_cnvkit_val
	
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
}

process freebayes {
	cpus 1
	time '40m'
	tag "$group"
	
	input:
		set group, id, type, file(bams), file(bais), file(bqsr) from bam_freebayes.groupTuple()
		each file(bed) from beds_freebayes

	output:
		set val("freebayes"), group, file("freebayes_${bed}.vcf") into vcfparts_freebayes

	when:
		params.freebayes

	script:
		dp = 500
		if (params.assay == "solid") {
			dp = 80
		}
			

		if( id.size() >= 2 ) {

			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 ${bams[tumor_idx]} ${bams[normal_idx]} > freebayes_${bed}.vcf.raw
			vcffilter -F LowCov -f "DP > $dp" -f "QA > 1500" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
			filter_freebayes_somatic.pl freebayes_${bed}.filt1.vcf ${id[tumor_idx]} ${id[normal_idx]} > freebayes_${bed}.vcf
			"""
		}
		else if( id.size() == 1 ) {
			"""
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bams > freebayes_${bed}.vcf.raw
			vcffilter -F LowCov -f "DP > $dp" -f "QA > 1500" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
			filter_freebayes_unpaired.pl freebayes_${bed}.filt1.vcf > freebayes_${bed}.vcf
			"""
		}
}


process vardict {
	cpus 1
	time '2h'
	tag "$group"
	memory '15GB'

	input:
		set group, id, type, file(bams), file(bais), file(bqsr) from bam_vardict.groupTuple()
		each file(bed) from beds_vardict

	output:
		set val("vardict"), group, file("vardict_${bed}.vcf") into vcfparts_vardict

	when:
		params.vardict
    
	script:
		if( id.size() >= 2 ) {

			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			vardict-java -G $genome_file -f 0.01 -N ${id[tumor_idx]} -b "${bams[tumor_idx]}|${bams[normal_idx]}" -c 1 -S 2 -E 3 -g 4 -U $bed \\
			| testsomatic.R | var2vcf_paired.pl -N "${id[tumor_idx]}|${id[normal_idx]}" -f 0.01 > vardict_${bed}.vcf.raw

			filter_vardict_somatic.pl vardict_${bed}.vcf.raw ${id[tumor_idx]} ${id[normal_idx]} > vardict_${bed}.vcf
			"""
		}
		else if( id.size() == 1 ) {
			"""
			vardict-java -G $genome_file -f 0.03 -N ${id[0]} -b ${bams[0]} -c 1 -S 2 -E 3 -g 4 -U $bed | teststrandbias.R | var2vcf_valid.pl -N ${id[0]} -E -f 0.01 > vardict_${bed}.vcf.raw
			filter_vardict_unpaired.pl vardict_${bed}.vcf.raw > vardict_${bed}.vcf
			"""
		}
}


process tnscope {
	cpus params.cpu_some
	time '2h'   
	tag "$group" 

	input:
		set group, id, type, file(bams), file(bais), file(bqsr) from bam_tnscope.groupTuple()
		each file(bed) from beds_tnscope

	output:
		set val("tnscope"), group, file("tnscope_${bed}.vcf") into vcfparts_tnscope

	when:
		params.tnscope

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

		if( id.size() >= 2 ) {
			"""
			sentieon driver -t ${task.cpus} \\
				-r $genome_file \\
				-i ${bams[tumor_idx]} -q ${bqsr[tumor_idx]} \\
				-i ${bams[normal_idx]} -q ${bqsr[normal_idx]} \\
				--interval $bed --algo TNscope \\
				--tumor_sample ${id[tumor_idx]} --normal_sample ${id[normal_idx]} \\
				--clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
				--min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
				tnscope_${bed}.vcf.raw

			filter_tnscope_somatic.pl tnscope_${bed}.vcf.raw ${id[tumor_idx]} ${id[normal_idx]} > tnscope_${bed}.vcf

			"""
		}
		else {
			"""
			sentieon driver -t ${task.cpus} -r $genome_file \\
				-i ${bams} -q ${bqsr} \\
				--interval $bed --algo TNscope \\
				--tumor_sample ${id[0]} \\
				--clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
				--min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
				tnscope_${bed}.vcf.raw

			filter_tnscope_unpaired.pl tnscope_${bed}.vcf.raw > tnscope_${bed}.vcf
			""" 
		}
}


process pindel {
	cpus params.cpu_some
	time '1h'
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	tag "$group"

	input:
		set group, id, type, file(bams), file(bais), file(ins_size) from all_pindel.groupTuple()

	output:
		set group, val("pindel"), file("${group}_pindel.vcf") into vcf_pindel

	when:
		params.pindel

	script:
		if( id.size() >= 2 ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			ins_tumor = ins_size[tumor_idx]
			ins_normal = ins_size[normal_idx]
			bam_tumor = bams[tumor_idx]
			bam_normal = bams[normal_idx]
			id_tumor = id[tumor_idx]
			id_normal = id[normal_idx]

			"""
			INS_T="\$(sed -n '3p' $ins_tumor | cut -f 1 | awk '{print int(\$1+0.5)}')"
			INS_N="\$(sed -n '3p' $ins_normal | cut -f 1 | awk '{print int(\$1+0.5)}')"
			echo "$bam_tumor\t\$INS_T\t$id_tumor" > pindel_config
			echo "$bam_normal\t\$INS_N\t$id_normal" >> pindel_config

			pindel -f $genome_file -w 0.1 -x 2 -i pindel_config -j $params.pindel_regions_bed -o tmpout -T ${task.cpus}
			pindel2vcf -P tmpout -r $genome_file -R hg19 -d 2015-01-01 -v ${group}_pindel_unfilt.vcf -is 10 -e 30 -he 0.01
			filter_pindel_somatic.pl ${group}_pindel_unfilt.vcf ${group}_pindel.vcf
			"""
		}
		else {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			ins_tumor = ins_size[tumor_idx]
			bam_tumor = bams[tumor_idx]
			id_tumor = id[tumor_idx]

			"""
			INS_T="\$(sed -n '3p' $ins_tumor | cut -f 1 | awk '{print int(\$1+0.5)}')"
			echo "$bam_tumor\t\$INS_T\t$id_tumor" > pindel_config

			pindel -f $genome_file -w 0.1 -x 2 -i pindel_config -j $params.pindel_regions_bed -o tmpout -T ${task.cpus}
			pindel2vcf -P tmpout -r $genome_file -R hg19 -d 2015-01-01 -v ${group}_pindel_unfilt.vcf -is 10 -e 30 -he 0.01
			filter_pindel_somatic.pl ${group}_pindel_unfilt.vcf ${group}_pindel.vcf
			"""
		}

}


// Prepare vcf parts for concatenation
vcfparts_freebayes = vcfparts_freebayes.groupTuple(by:[0,1])
vcfparts_tnscope   = vcfparts_tnscope.groupTuple(by:[0,1])
vcfparts_vardict   = vcfparts_vardict.groupTuple(by:[0,1])
vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict).mix(vcfparts_tnscope)

process concatenate_vcfs {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	time '20m'    
	tag "$group"

	input:
		set vc, group, file(vcfs) from vcfs_to_concat

	output:
		set group, vc, file("${group}_${vc}.vcf.gz") into concatenated_vcfs, vcf_cnvkit

	"""
	vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt normalize ${vc}.decomposed.vcf.gz -r $genome_file | vt uniq - -o ${group}_${vc}.vcf.gz
	"""
}


process cnvkit {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true, pattern: '*.vcf'
	publishDir "${OUTDIR}/gens", mode: 'copy', overwrite: true, pattern: '*.bed.gz*'
	publishDir "${CRONDIR}/gens", mode: 'copy', overwrite: true, pattern: '*.gens'
	publishDir "${OUTDIR}/cnvkit", mode: 'copy', overwrite: true, pattern: '*.cnvkit'
	publishDir "${OUTDIR}/cnvkit/HRD", mode: 'copy', overwrite: true, pattern: '*.HRD'
	cpus 1
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = '/fs1/resources/containers/cnvkit099.sif'
	
	input:
		set gr, id, type, file(bam), file(bai), file(bqsr), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV), vc, file(vcf) from bam_cnvkit.join(qc_cnvkit_val, by:[0,1]) \
			.combine(vcf_cnvkit.filter { item -> item[1] == 'freebayes' }, by:[0])
		
	output:
		set gr, id, type, file("${gr}.${id}.cnvkit_overview.png"), file("${gr}.${id}.call.cns"), file("${gr}.${id}.cnr"), file("${gr}.${id}.filtered") into geneplot_cnvkit
		set gr, id, type, file("${gr}.${id}.filtered.vcf") into cnvkit_vcf 
		file("${gr}.${id}.cns") into cns_notcalled
		file("*.bed.gz*")
		file("${id}.gens") into gens_middleman
		file("${gr}.${id}_logr_ballele.cnvkit")
		file("${gr}.${id}_seg.cnvkit")
		file("${id}.HRD")

	when:
		params.cnvkit

	script:
		freebayes_idx = vc.findIndexOf{ it == 'freebayes' }

	"""
	set +eu
	source activate py2
	set -eu

	cnvkit.py batch $bam -r $params.cnvkit_reference -d results/
	cnvkit.py call results/*sort.cns -v $vcf -o ${gr}.${id}.call.cns
	filter_cnvkit.pl ${gr}.${id}.call.cns $MEAN_DEPTH 1000000 > ${gr}.${id}.filtered
	cnvkit.py export vcf ${gr}.${id}.filtered -i "$id" > ${gr}.${id}.filtered.vcf		
	cnvkit.py scatter -s results/*sort.cn{s,r} -o ${gr}.${id}.cnvkit_overview.png -v ${vcf[freebayes_idx]} -i $id
	cp results/*sort.cnr ${gr}.${id}.cnr
	cp results/*sort.cns ${gr}.${id}.cns
	cnvkit.py export nexus-ogt -o 11-FF-HRD-GMSSTv1-0-210203.nexus -o ${gr}.${id}_logr_ballele.cnvkit ${gr}.${id}.cnr ${vcf[freebayes_idx]}
	cnvkit.py export seg -o ${gr}.${id}_seg.cnvkit ${gr}.${id}.cns
	generate_gens_data_from_cnvkit.pl ${gr}.${id}.cnr $vcf $id
	python /fs1/viktor/SomaticPanelPipeline/bin/HRD.py ${gr}.${id}.call.cns tmp > ${id}.HRD
	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz" > ${id}.gens
	"""
}

// Plot specific gene-regions. CNVkit 0.9.6 and forward introduced a bug in region plot, use 0.9.5 (091 wrong name, container has 0.9.5)
process gene_plot {
	publishDir "${OUTDIR}/plots", mode: 'copy', overwrite: true, pattern: '*.png'
	cpus 1
	time '5m'
	tag "$id"

	input:
		set gr, id, type, file(overview), file(cns), file(cnr), file(filtered) from geneplot_cnvkit

	output:
		set gr, id, type, file("${gr}.${id}.cnvkit.png") into cnvplot_coyote

	script:

		if (params.assay == "PARP_inhib") {
			"""
            set +eu
            source activate old-cnvkit
            set -eu
			cnvkit.py scatter -s $cns $cnr -c 13:32165479-32549672 -o brca2.png --title 'BRCA2'
			cnvkit.py scatter -s $cns $cnr -c 17:42894294-43350132 -o brca1.png --title 'BRCA1'
			montage -mode concatenate -tile 1x *.png ${gr}.${id}.cnvkit.png
			"""		
		}
		else {
			"""
			mv ${gr}.${id}.cnvkit_overview.png ${gr}.${id}.cnvkit.png
			"""
		}


}

process melt {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 2
	//container = '/fs1/resources/containers/container_twist-brca.sif'
	memory '50 GB'
	tag "$group"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'


	input:
		set group, id, type, file(bam), file(bai), file(bqsr), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_melt \
			.join(qc_melt_val, by: [0,1,2])
		
	when:
		params.melt

	output:
		set group, file("${id}_${type}.melt.merged.vcf") into melt_vcf

	"""
	set +eu
	source activate java8
	set -eu
	java -jar  /opt/MELT.jar Single \\
		-bamfile $bam \\
		-r 150 \\
		-h $genome_file \\
		-n $params.bed_melt \\
		-z 50000 \\
		-d 50 -t $params.mei_list \\
		-w . \\
		-b 1/2/3/4/5/6/7/8/9/10/11/12/14/15/16/18/19/20/21/22 \\
		-c $MEAN_DEPTH \\
		-cov $COV_DEV \\
		-e $INS_SIZE
        source deactivate
	merge_melt.pl $params.meltheader $id
	mv ${id}.melt.merged.vcf ${id}_${type}.melt.merged.vcf
	"""

}

// MANTA SINGLE AND PAIRED
process manta {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 16
	time '10h'
	//container = '/fs1/resources/containers/wgs_2020-03-25.sif'
	tag "$group"
	scratch true
	memory '10GB'
	stageInMode 'copy'
	stageOutMode 'copy'
	
	input:
		set group, id, type, file(bam), file(bai), file(bqsr) from bam_manta.groupTuple()

	output:
		set group, file("${group}_manta.vcf") into manta_vcf

	when:
		params.manta
	
	script:
		if(id.size() >= 2) { 
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			normal = bam[normal_idx]
			normal_id = id[normal_idx]
			tumor = bam[tumor_idx]
			tumor_id = id[tumor_idx]

			"""
            set +eu
            source activate py2
            set -eu
			configManta.py \\
				--tumorBam $tumor \\
				--normalBam $normal \\
				--reference $genome_file \\
				--exome \\
				--callRegions $params.bedgz \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			#filter_manta_paired.pl results/variants/somaticSV.vcf.gz > ${group}_manta.vcf
			mv results/variants/somaticSV.vcf.gz ${group}_manta.vcf.gz
			gunzip ${group}_manta.vcf.gz
			"""
		}
		else {
			"""
            set +eu
            source activate py2
            set -eu
			configManta.py \\
				--tumorBam $bam \\
				--reference $genome_file \\
				--exome \\
				--callRegions $params.bedgz \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			#filter_manta.pl results/variants/tumorSV.vcf.gz > ${group}_manta.vcf
			mv results/variants/tumorSV.vcf.gz ${group}_manta.vcf.gz
			gunzip ${group}_manta.vcf.gz
			"""
		}
}

// Delly SINGLE AND PAIRED
process delly {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 2
	time '20h'
	memory '10GB'
	//container = '/fs1/resources/containers/wgs_2020-03-25.sif'
	tag "$group"
		
	input:
		set group, id, type, file(bam), file(bai), file(bqsr) from bam_delly.groupTuple()

	output:
		set group, file("${group}.delly.filtered.vcf") into delly_vcf

	when:
		params.delly
	
	script:
		if(id.size() >= 2) { 
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			normal = bam[normal_idx]
			normal_id = id[normal_idx]
			tumor = bam[tumor_idx]
			tumor_id = id[tumor_idx]

			"""
			delly call -g $genome_file -o ${group}.delly.bcf $tumor $normal
			bcftools view ${group}.delly.bcf > ${group}.delly.vcf
			filter_delly.pl --vcf ${group}.delly.vcf --bed $params.regions_bed > ${group}.delly.filtered.vcf
			"""
		}
		else {
			"""
			delly call -g $genome_file -o ${group}.delly.bcf $bam
			bcftools view ${group}.delly.bcf > ${group}.delly.vcf
			filter_delly.pl --vcf ${group}.delly.vcf --bed $params.regions_bed > ${group}.delly.filtered.vcf
			"""
		}
}

process aggregate_CNVkit {
    time '20m'
    tag "$group"
	cpus 2
	memory '1GB'
	//publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true

    input:
        set group, id, type, file(vcf) from cnvkit_vcf.groupTuple()
       
    output:
        set group, file("${group}_cnvkit_agg.vcf") into cnvkit_vcfagg
       
    script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

	if (id.size() >= 2 ) {
		"""
		aggregate_CNVkit.pl ${vcf[tumor_idx]} ${id[tumor_idx]} ${vcf[normal_idx]} ${id[normal_idx]} > ${group}_cnvkit_agg.vcf
		"""
	}
	else {
		"""
		mv $vcf ${group}_cnvkit_agg.vcf
		"""
	}

}

process gatk_coverage {
    cpus 10
    memory '50GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    tag "$id"   

	when:
		params.gatk_cnv

    input:
        set group, id, type, file(bam), file(bai) from bam_gatk

    output:
        set group, id, file("${id}.tsv") into call_ploidy, call_cnv

    """
	export THEANO_FLAGS="base_compiledir=."
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
	set +u
	source activate gatk
    gatk --java-options "-Xmx20g" CollectReadCounts \\
        -L $params.gatk_intervals \\
        -R $params.genome_file \\
        -imr OVERLAPPING_ONLY \\
        -I $bam \\
        --format TSV -O ${id}.tsv
    """
}

process gatk_call_ploidy {
    cpus 10
    memory '40GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    tag "$id"

    input:
        set group, id, file(tsv) from call_ploidy

    output:
        set group, id, file("ploidy.tar") into ploidy_to_cnvcall, ploidy_to_post

    """
	export THEANO_FLAGS="base_compiledir=."
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
	set +u
	source activate gatk
    gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
        --model $params.ploidymodel \\
        -I $tsv \\
        -O ploidy/ \\
        --output-prefix $group
    tar -cvf ploidy.tar ploidy/
    """
}

process gatk_call_cnv {
    cpus 8
    memory '45GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    tag "$id"

    input:
        set group, id, file(tsv), file(ploidy), i, refpart \
            from call_cnv.join(ploidy_to_cnvcall, by: [0,1]).combine(gatk_ref)

    output:
        set group, id, i, file("${group}_${i}.tar") into postprocessgatk

    """
	export THEANO_FLAGS="base_compiledir=."
	set +u
	source activate gatk
	export HOME=/local/scratch
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    tar -xvf ploidy.tar
    mkdir ${group}_${i}
    gatk --java-options "-Xmx25g" GermlineCNVCaller \\
        --run-mode CASE \\
        -I $tsv \\
        --contig-ploidy-calls ploidy/${group}-calls/ \\
        --model ${refpart} \\
        --output ${group}_${i}/ \\
        --output-prefix ${group}_${i}
    tar -cvf ${group}_${i}.tar ${group}_${i}/
    """
}

process postprocessgatk {
    cpus 8
    memory '40GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	
    //scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
    publishDir "${OUTDIR}/gatk_cnv/", mode: 'copy', overwrite: 'true'
    tag "$id"


    input:
        set group, id, i, file(tar), file(ploidy), shard_no, shard \
			from postprocessgatk.groupTuple(by: [0,1]).join(ploidy_to_post, by: [0,1]).combine(gatk_postprocess.groupTuple(by: [3]))


    output:
        set group, id, \
            file("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), \
            file("genotyped-segments-${group}-vs-cohort30.vcf.gz"), \
            file("denoised-${group}-vs-cohort30.vcf.gz") into called_gatk

    script:
        modelshards = shard.join(' --model-shard-path ') // join each reference shard
        caseshards = []
        for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
            tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls' 
            caseshards = caseshards + tmp
        }
        caseshards = caseshards.join( ' --calls-shard-path ')
 	shell:
	'''
	THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
	for model in !{tar}; do
	tar -xvf $model
	done
    tar -xvf !{ploidy}
	set +u
	source activate gatk
    export MKL_NUM_THREADS=!{task.cpus}
    export OMP_NUM_THREADS=!{task.cpus}
    gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
        --allosomal-contig X --allosomal-contig Y \
        --contig-ploidy-calls ploidy/!{group}-calls/ \
        --sample-index 0 \\
        --output-genotyped-intervals genotyped-intervals-!{group}-vs-cohort30.vcf.gz \
        --output-genotyped-segments genotyped-segments-!{group}-vs-cohort30.vcf.gz \
        --output-denoised-copy-ratios denoised-!{group}-vs-cohort30.vcf.gz \
        --sequence-dictionary !{params.GENOMEDICT} \
        --calls-shard-path !{caseshards} \
        --model-shard-path !{modelshards}
	'''

}

process filter_merge_gatk {
	cpus 1
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "${OUTDIR}/gatk_cnv", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(inter), file(gatk), file(denoised) from called_gatk

	output:
		set group, id, file("${id}.gatk.filtred.merged.vcf") into merged_gatk

	"""
	filter_gatk.pl $gatk > ${id}.gatk.filtered.vcf
	mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtred.merged.vcf
	"""
}



process aggregate_cnvs {
        cpus 2
        memory '1GB'
        publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
        time '20m'
        tag "$group"

		when:
			!params.single_cnvcaller

        input:
			set group, file(vcfs) from cnvkit_vcfagg.mix(manta_vcf,delly_vcf).groupTuple().view()
			
        output:
            set group, file("${group}.merged.vcf") into cnvs

		script:
			if (vcfs.size() > 1) {
				// for each sv-caller add idx, find vcf and find priority, add in priority order! //
				// index of vcfs added from mix //
				manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
				delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
				cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }

				// find vcfs //
				manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
				delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
				cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
				tmp = manta + delly + cnvkit
				tmp = tmp - null
				vcfs_svdb = tmp.join(' ')

				// find priorities //
				mantap = manta_idx >= 0 ? 'manta' : null
				dellyp = delly_idx >= 0 ? 'delly' : null
				cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
				tmpp = [mantap, dellyp, cnvkitp]
				tmpp = tmpp - null
				priority = tmpp.join(',')

				"""
				svdb --merge --vcf $vcfs_svdb --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority $priority > ${group}.merged.vcf
				"""

			}
			else {
				"""
				cat $vcfs > ${group}.merged.vcf
				"""
			}

}

process standardize_cnvs {
	cpus 1
	memory '1GB'
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	time '20m'
	tag "$group"

	input:
		set group, vcf, id, type, tissue from cnvs.mix(melt_vcf.groupTuple()).groupTuple().join(meta_cnvkit.groupTuple()).view()
		
	
	output:
		set group, file("${group}.cnvs.agg.vcf") into standard_cnvs
	
	script:
		// check vcf.size, if two elements then melt was run, vcf = agg,melt(tumor) //
		if (vcf.size() > 1) {
			cnv_idx = vcf.findIndexOf{ it =~ 'agg' }
			melt_idx = cnv_idx == 0 ? 1 : 0
			melt_t = vcf[melt_idx].findIndexOf{ it =~ 'T.melt' }
			melt = vcf[melt_idx][melt_t]
			tmp = [vcf[cnv_idx], melt]
			vcf = tmp.join(',')
		}
		else {
			vcf = vcf.join(',')
		}
		// if paired
		if( id.size() >= 2 ) {
				tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
				normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

				"""
				aggregate_cnv2_vcf.pl --vcfs $vcf \\
						--tumor-id ${id[tumor_idx]} \\
						--normal-id ${id[normal_idx]} \\
						--paired paired \\
						--sample-order ${id[tumor_idx]},${id[normal_idx]} > ${group}.cnvs.agg.vcf
				"""
		}
		// tumor only
		else {

			"""
			aggregate_cnv2_vcf.pl --vcfs $vcf --paired no > ${group}.cnvs.agg.vcf
			"""
		}
}

process single_cnv_pipe {
       time '2m'
       tag "$group"

       when:
               params.single_cnvcaller

       input:
               set group, id, type, file(read1), file(read2) from meta_nocnv
       
       output:
               set group, file("${group}.cnvs.agg.vcf") into cnvs_singlecaller
       
       script:
       """
       echo singe_cnv_caller_pipeline > ${group}.cnvs.agg.vcf
       """
}

process aggregate_vcfs {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	time '20m'
	tag "$group"

	input:
		set group, vc, file(vcfs), id, type, tissue, file(cnvs) from concatenated_vcfs.mix(vcf_pindel).groupTuple().join(meta_aggregate.groupTuple()).join(standard_cnvs.mix(cnvs_singlecaller))

	output:
		set group, file("${group}.agg.vcf") into vcf_pon, vcf_done

	script:
		sample_order = id[0]
		if( id.size() >= 2 ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			sample_order = id[tumor_idx]+","+id[normal_idx]
		}
		if (params.single_cnvcaller) {
			"""
			aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.vcf
			"""
		}
		else {
			"""
			aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.tmp.vcf
			vcf-concat ${group}.agg.tmp.vcf $cnvs | vcf-sort -c > ${group}.agg.vcf
			"""
		}

}

process pon_filter {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 1
	time '1h'
	tag "$group"
	memory '32 GB'

	input:
		set group, file(vcf), id, type, tissue from vcf_pon.join(meta_pon.groupTuple())
		
	output:
		set group, file("${group}.agg.pon.vcf") into vcf_vep

	script:
	if (params.assay == 'myeloid') {
			def pons = []
			if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
			if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
			if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
			def pons_str = pons.join(",")
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }

		"""
		filter_with_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${id[tumor_idx]} > ${group}.agg.pon.vcf
		"""
	}
	// weird placement, no PON for PARP_inhib, Adds enigma-db to vcf. Move to separate process?
	else if (params.assay == 'PARP_inhib') {
		"""
		vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno $vcf > ${group}.agg.pon.vcf
		"""
	}
	else {
		"""
		mv $vcf ${group}.agg.pon.vcf
		"""
	}
}

process annotate_vep {
	container = params.vepcon
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus params.cpu_many
	time '1h'
	tag "$group"
    
	input:
		set group, file(vcf) from vcf_vep
    
	output:
		set group, file("${group}.agg.pon.vep.vcf") into vcf_germline

	"""
	vep -i ${vcf} -o ${group}.agg.pon.vep.vcf \\
	--offline --merged --everything --vcf --no_stats \\
	--fork ${task.cpus} \\
	--force_overwrite \\
	--plugin CADD $params.CADD --plugin LoFtool \\
	--fasta $params.VEP_FASTA \\
	--dir_cache $params.VEP_CACHE --dir_plugins $params.VEP_CACHE/Plugins \\
	--distance 200 \\
	--custom $params.GNOMAD,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
	--custom $params.COSMIC,COSMIC,vcf,exact,0,CNT \\
	--cache \\
	"""
}

process mark_germlines {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus params.cpu_many
	time '20m'
	tag "$group"

	input:
		set group, file(vcf), id, type, tissue from vcf_germline.join(meta_germline.groupTuple())

		
	output:
		set group, file("${group}.agg.pon.vep.markgerm.vcf") into vcf_umi


	script:
		if( id.size() >= 2 ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			mark_germlines.pl --vcf $vcf --tumor-id ${id[tumor_idx]} --normal-id ${id[normal_idx]} --assay $params.assay > ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
		else if( id.size() == 1 ) {
			"""
			mark_germlines.pl --vcf $vcf --tumor-id ${id[0]} --assay $params.assay > ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
}

process umi_confirm {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 2
	time '8h'
	tag "$group"

	when:
		params.umi

	input:
		set group, file(vcf), id, type, file(bam), file(bai) from vcf_umi.join(bam_umi_confirm.groupTuple())
	
	output:
		set group, file("${group}.agg.pon.vep.markgerm.umi*") into vcf_coyote


	script:
		if (params.conform) {
	
			if( id.size() >= 2 ) {
				tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
				normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }

				"""
				UMIconfirm_vcf.py ${bam[tumor_idx]} $vcf $genome_file ${id[tumor_idx]} > umitmp.vcf
				UMIconfirm_vcf.py ${bam[normal_idx]} umitmp.vcf $genome_file ${id[normal_idx]} > ${group}.agg.pon.vep.markgerm.umi.vcf
				"""
			}
			else if( id.size() == 1 ) {
				tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }

				"""
				UMIconfirm_vcf.py ${bam[tumor_idx]} $vcf $genome_file ${id[tumor_idx]} > ${group}.agg.pon.vep.markgerm.umi.vcf
				"""
			}

		}
		else {
			"""
			cp $vcf ${group}.agg.pon.vep.markgerm.umino.vcf
			"""
		}
}


process coyote {
	publishDir "${params.crondir}/coyote", mode: 'copy', overwrite: true
	cpus 1
	time '10m'
	tag "$group"

	input:
		set group, file(vcf),  type, lims_id, pool_id, id, cnv_type, \
			file(cnvplot), tissue_c, lowcov_type, file(lowcov) from \
			vcf_coyote.join(meta_coyote.groupTuple()).join(cnvplot_coyote.join(meta_cnvplot, by:[0,1,2]).groupTuple()).join(lowcov_coyote.groupTuple())


	output:
		file("${group}.coyote")

	when:
		!params.noupload

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
		tumor_idx_cnv = cnv_type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx_cnv = cnv_type.findIndexOf{ it == 'normal' || it == 'N' }
		cnv_index = tumor_idx_cnv
		tumor_idx_lowcov = lowcov_type.findIndexOf{ it == 'tumor' || it == 'T' }
		if( id.size() >= 2 ) {
			group = group + 'p'
		}


	"""
	echo "import_myeloid_to_coyote_vep_gms.pl --group $params.coyote_group \\
		--vcf /access/${params.subdir}/vcf/${vcf} --id ${group} \\
		--cnv /access/${params.subdir}/plots/${cnvplot[cnv_index]} \\
		--clarity-sample-id ${lims_id[tumor_idx]} \\
		--lowcov /access/${params.subdir}/QC/${lowcov[tumor_idx_lowcov]} \\
                --build 38 \\
                --gens ${group} \\
		--clarity-pool-id ${pool_id[tumor_idx]}" > ${group}.coyote
	"""
}
