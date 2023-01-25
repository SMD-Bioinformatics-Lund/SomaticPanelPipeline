process MANTA {
	publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: true, pattern: '*.vcf'
	cpus 16
	time '10h'
	tag "$group"
	//scratch true
	memory '10GB'
	//stageInMode 'copy'
//	stageOutMode 'copy'
	
	input:
		tuple val(group), val(meta), file(bam), file(bai)

	output:
		tuple val(group), file("${meta.id[tumor_idx]}_manta.vcf"), emit: manta_vcf_tumor
		tuple val(group), file("${meta.id[tumor_idx]}_bnd_manta.vcf"), emit: manta_vcf_bnd_tumor
		tuple val(group), file("${meta.id[normal_idx]}_manta.vcf"), optional: true, emit: manta_vcf_normal
		tuple val(group), file("${meta.id[normal_idx]}_bnd_manta.vcf"), optional: true, emit: manta_vcf_bnd_normal
		

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
				--callRegions $params.bedgz \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			gunzip results/variants/somaticSV.vcf.gz
			gunzip results/variants/diploidSV.vcf.gz
			grep -v BND results/variants/somaticSV.vcf > ${meta.id[tumor_idx]}_manta.vcf
			grep ^# results/variants/somaticSV.vcf > ${meta.id[tumor_idx]}_bnd_manta.vcf
			grep -v ^# results/variants/somaticSV.vcf | grep BND >> ${meta.id[tumor_idx]}_bnd_manta.vcf
			grep -v BND results/variants/diploidSV.vcf > ${meta.id[normal_idx]}_manta.vcf
			grep ^# results/variants/diploidSV.vcf > ${meta.id[normal_idx]}_bnd_manta.vcf
			grep -v ^# results/variants/diploidSV.vcf | grep BND >> ${meta.id[normal_idx]}_bnd_manta.vcf
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
				--callRegions $params.bedgz \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			gunzip results/variants/tumorSV.vcf.gz
			grep -v BND results/variants/tumorSV.vcf > ${meta.id[tumor_idx]}_manta.vcf
			grep ^# results/variants/tumorSV.vcf > ${meta.id[tumor_idx]}_bnd_manta.vcf
			grep -v ^# results/variants/tumorSV.vcf | grep BND >> ${meta.id[tumor_idx]}_bnd_manta.vcf
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
			touch ${meta.id[tumor_idx]}_manta.vcf
			touch ${meta.id[normal_idx]}_manta.vcf
			"""
		}
		else {
			"""
			touch ${meta.id[tumor_idx]}_manta.vcf
			"""
		}

}