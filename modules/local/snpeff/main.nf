process SNPEFF {
 	publishDir "${params.outdir}/${params.subdir}/fusions/", mode: 'copy', overwrite: true, pattern: '*.vcf'
	cpus 16
	time '10h'
	tag "$group"
	// scratch true
	memory '4GB'
	// stageInMode 'copy'
	// stageOutMode 'copy'
    container = "/fs1/resources/containers/snpeff_4.3.1t2.sif"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("${group}.merged.annotated.vcf"), emit: snpeff_vcf

    script:
        """
        snpEff -Xmx4g -configOption data.dir=/fs1/resources/ref/hg38/snpeff/ GRCh38.86 ${vcf} > ${group}.merged.annotated.vcf
        """
    stub:
        """
        touch ${group}.merged.annotated.vcf
        """ 
}