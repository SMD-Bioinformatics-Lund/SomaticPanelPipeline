process CONCATENATE_VCFS {
	cpus 1
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	time '20m'    
	tag "$group"

	input:
		tuple val(vc), val(group), file(vcfs)

	output:
		tuple val(group), val(vc), file("${group}_${vc}.vcf.gz"), emit: concatenated_vcfs

	"""
	vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt normalize ${vc}.decomposed.vcf.gz -r $params.genome_file | vt uniq - -o ${group}_${vc}.vcf.gz
	"""
}

process AGGREGATE_VCFS {
	cpus 1
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	time '20m'
	tag "$group"

	input:
		tuple val(group), val(vc), file(vcfs), val(id), val(type), val(FFPE)
		//row.group, row.id, row.type, (row.containsKey("ffpe") ? row.ffpe

	output:
		tuple val(group), file("${group}.agg.vcf"), emit: vcf_concat

	script:
		sample_order = id[0]
		if( id.size() >= 2 ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			println(id[tumor_idx])
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

process CONTAMINATION {
	publishDir "${params.outdir}/${params.subdir}/QC/contamination", mode: 'copy', overwrite: true, pattern: "*.png"
	publishDir "${params.outdir}/${params.subdir}/QC/contamination", mode: 'copy', overwrite: true, pattern: "*.txt"
	publishDir "${params.crondir}/contamination", mode: 'copy', overwrite: true, pattern: "*.contamination"
	container = "/fs1/resources/containers/perl-gd.sif"
	//errorStrategy 'ignore'
	cpus 1
	time '10m'
	tag "$group"

	input:
		tuple val(group), file(vcf), val(id), val(type), val(sequencing_run)

	output:
		tuple val(group), file("*.txt"), file("*.png"), emit: contamination_result_files
		tuple val(group), file("*.contamination"), emit: contamination_cdm

	script:
		if(id.size() >= 2) { 
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			find_contaminant.pl --vcf $vcf --case-id ${id[tumor_idx]} --assay ${params.cdm} --detect-level 0.01 > ${id[tumor_idx]}.value
			echo "--overwrite --sample-id ${id[tumor_idx]} --run-folder ${sequencing_run[tumor_idx]} --assay ${params.cdm} --contamination" > ${tumor_id}.1
			paste -d " " ${id[tumor_idx]}.1 ${id[tumor_idx]}.value > ${id[tumor_idx]}.contamination
			find_contaminant.pl --vcf $vcf --case-id ${id[tumor_idx]} --assay ${params.cdm} --detect-level 0.01 --normal > ${id[normal_idx]}.value
			echo "--overwrite --sample-id ${id[normal_idx]} --run-folder ${sequencing_run[normal_idx]} --assay ${params.cdm} --contamination" > ${id[normal_idx]}.1
			paste -d " " ${id[normal_idx]}.1 ${id[normal_idx]}.value > ${id[normal_idx]}.contamination
			"""
		}
		else {
			"""
			find_contaminant.pl --vcf $vcf --case-id ${id[0]} --assay ${params.cdm} --detect-level 0.01 > ${id[0]}.value
			echo "--overwrite --sample-id ${id[0]} --run-folder ${sequencing_run[0]} --assay ${params.cdm} --contamination" > ${id[0]}.1
			paste -d " " ${id[0]}.1 ${id[0]}.value > ${id[0]}.contamination
			"""
		}


}