process FREEBAYES {
	cpus 1
	time '40m'
	tag "$group"
	
	input:
		tuple val(group), val(meta), file(bams), file(bais), file(bqsr)
		each file(bed)

	output:
		tuple val("freebayes"),  val(group), file("freebayes_${bed}.vcf"), emit: vcfparts_freebayes

	when:
		params.freebayes

	script:
		dp = 500
		if (params.assay == "solid") {
			dp = 80
		}

		if( meta.id.size() >= 2 ) {

			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

			"""
			freebayes -f $params.genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 ${bams[tumor_idx]} ${bams[normal_idx]} > freebayes_${bed}.vcf.raw
			vcffilter -F LowCov -f "DP > $dp" -f "QA > 1500" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
			filter_freebayes_somatic.pl freebayes_${bed}.filt1.vcf ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > freebayes_${bed}.vcf
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			freebayes -f $params.genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bams > freebayes_${bed}.vcf.raw
			vcffilter -F LowCov -f "DP > $dp" -f "QA > 1500" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
			filter_freebayes_unpaired.pl freebayes_${bed}.filt1.vcf > freebayes_${bed}.vcf
			"""
		}
	stub:
	if( meta.id.size() >= 2 ) {
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
		"""
		echo tumor:${bams[tumor_idx]} ${meta.id[tumor_idx]} normal:${bams[normal_idx]} ${meta.id[normal_idx]}
		touch freebayes_${bed}.vcf
		"""
	}
	else {
		"""
		echo tumor:$bams
		touch freebayes_${bed}.vcf
		"""
	}

}