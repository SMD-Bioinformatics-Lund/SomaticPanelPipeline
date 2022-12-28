process MANTA {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 16
	time '10h'
	tag "$group"
//	scratch true
	memory '10GB'
	// stageInMode 'copy'
	// stageOutMode 'copy'
	
	input:
		tuple val(group), val(meta), file(bam), file(bai)

	output:
		tuple val(group), file("${group}_manta.vcf"), emit: manta_vcf

	when:
		params.manta
	
	script:
		if(meta.id.size() == 2) { 
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			normal = bam[normal_idx]
			normal_id = meta.id[normal_idx]
			tumor = bam[tumor_idx]
			tumor_id = meta.id[tumor_idx]

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
				--reference ${params.genome_file} \\
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