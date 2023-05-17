process MANTA {
	publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: true, pattern: '*.vcf'
	cpus 16
	time '10h'
	tag "$group"
	scratch true
	memory '10GB'
	stageInMode 'copy'
	stageOutMode 'copy'
	
	input:
		tuple val(group), val(meta), file(bam), file(bai)
		val(reference)
		val(type)

	output:
		tuple val(group), file("${meta.id[tumor_idx]}_manta.${type}.vcf"), emit: manta_vcf_tumor
		tuple val(group), file("${meta.id[normal_idx]}_manta.${type}.vcf"), optional: true, emit: manta_vcf_normal
		

	when:
		params.manta
	
	script:
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
		normal = bam[normal_idx]
		normal_id = meta.id[normal_idx]
		tumor = bam[tumor_idx]
		tumor_id = meta.id[tumor_idx]
		if(meta.id.size() == 2) { 
			"""
            set +eu
            source activate py2
            set -eu
			configManta.py \\
				--tumorBam $tumor \\
				--normalBam $normal \\
				--reference ${params.genome_file} \\
				--exome \\
				--callRegions $reference \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			mv results/variants/somaticSV.vcf.gz ${meta.id[tumor_idx]}_manta.${type}.vcf.gz
			mv results/variants/diploidSV.vcf.gz ${meta.id[normal_idx]}_manta.${type}.vcf.gz
			gunzip ${meta.id[tumor_idx]}_manta.${type}.vcf.gz
			gunzip ${meta.id[normal_idx]}_manta.${type}.vcf.gz
			"""
		}
		else {
			"""
            set +eu
            source activate py2
            set -eu
			configManta.py \\
				--tumorBam $bam \\
				--reference ${params.genome_file} \\
				--exome \\
				--callRegions $reference \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			mv results/variants/tumorSV.vcf.gz ${meta.id[tumor_idx]}_manta.${type}.vcf.gz
			gunzip ${meta.id[tumor_idx]}_manta.${type}.vcf.gz
			"""
		}
	stub:
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
		normal = bam[normal_idx]
		normal_id = meta.id[normal_idx]
		tumor = bam[tumor_idx]
		tumor_id = meta.id[tumor_idx]
		if(meta.id.size() == 2) {
			"""
			touch ${meta.id[tumor_idx]}_manta.${type}.vcf ${meta.id[normal_idx]}_manta.${type}.vcf
			"""
		}
		else {
			"""
			touch ${meta.id[tumor_idx]}_manta.${type}.vcf ${meta.id[tumor_idx]}_manta_filtered.${type}.vcf
			"""
		}

}