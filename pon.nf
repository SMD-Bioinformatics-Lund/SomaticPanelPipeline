
params.csv = "/fs1/ram/Testing/SomaticPanel_LOD/PON/pon_creation_samples.csv"
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
	logFile = file("/fs1/ram/Testing/SomaticPanel_LOD/PON/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}

// Print active container
container = file(params.container).toRealPath()
println("container: "+container)

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.bam), file(row.bai)) }
    .set { sorted_umi_bams }


process bqsr_umi {
	cpus params.cpu_some
	memory '16 GB'
	time '1h'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set id, file(bam), file(bai) from sorted_umi_bams

	output:
		set id, file(bam), file(bai), file("${id}.bqsr.table") into bam_freebayes, bam_vardict, bam_tnscope
	when:
		params.umi

	"""
	sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo QualCal ${id}.bqsr.table
	"""
}

process freebayes {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 15
	time '2h'
	tag "$id"
	
	input:
		set id, file(bam), file(bai), file(bqsr) from bam_freebayes

	output:
		set id, file("${id}.freebayes.vcf.gz") into vcf_freebayes

	when:
		params.freebayes

	script:
        """
		freebayes -f $genome_file -t $params.regions_bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F $params.freebayes_var_freq_cutoff_up $bam > $id".freebayes.vcf.raw"
		vcffilter -F LowCov -f "DP > 500" -f "QA > 1500" $id".freebayes.vcf.raw" | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > $id".freebayes.filt1.vcf"
		filter_freebayes_unpaired.pl $id".freebayes.filt1.vcf" | vcf-sort -c | gzip -c > $id".freebayes.filtered.vcf.gz"

		vt decompose $id".freebayes.filtered.vcf.gz" -o $id".freebayes.decomposed.vcf.gz"
		vt normalize $id".freebayes.decomposed.vcf.gz" -r $genome_file | vt uniq - -o $id".freebayes.vcf.gz"
		"""
}

process vardict {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 15
	time '2h'
	tag "$id"
	memory '45GB'

	input:
		set id, file(bam), file(bai), file(bqsr) from bam_vardict

	output:
		set id, file("${id}.vardict.vcf.gz") into vcf_vardict

	when:
		params.vardict
    
	script:

		"""
		vardict-java -G $genome_file -f $params.vardict_var_freq_cutoff_up -N $id -b $bam -c 1 -S 2 -E 3 -g 4 -U $params.regions_bed | teststrandbias.R | var2vcf_valid.pl -N $id -E -f 0.01 > $id".vardict.vcf.raw"
		filter_vardict_unpaired.pl $id".vardict.vcf.raw" | vcf-sort -c | gzip -c > $id".vardict.filtered.vcf.gz"
	
		vt decompose $id".vardict.filtered.vcf.gz" -o $id".vardict.decomposed.vcf.gz"
		vt normalize $id".vardict.decomposed.vcf.gz" -r $genome_file | vt uniq - -o $id".vardict.vcf.gz"
		"""
}


process tnscope {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 15
	time '2h'   
	tag "$id" 

	input:
		set id, file(bam), file(bai), file(bqsr) from bam_tnscope

	output:
		set id, file("${id}.tnscope.vcf.gz") into vcf_tnscope

	when:
		params.tnscope

	script:

		"""
		sentieon driver -t ${task.cpus} -r $genome_file \\
			-i $bam -q $bqsr \\
			--interval $params.regions_bed --algo TNscope \\
			--tumor_sample $id \\
			--clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 \\
			--min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.0005 \\
			$id".tnscope.vcf.raw"
		filter_tnscope_unpaired.pl $id".tnscope.vcf.raw" | vcf-sort -c | gzip -c > $id".tnscope.filtered.vcf.gz"
	
		vt decompose $id".tnscope.filtered.vcf.gz" -o $id".tnscope.decomposed.vcf.gz"
		vt normalize $id".tnscope.decomposed.vcf.gz" -r $genome_file | vt uniq - -o $id".tnscope.vcf.gz"
		""" 

}