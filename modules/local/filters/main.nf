process PON_FILTER {
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(meta), file(vcf) 
		
	output:
		tuple val(group), val(meta), file("${group}.agg.pon.vcf"), emit: vcf_pon

	script:
        def pons = []
        if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
        if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
        if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
        def pons_str = pons.join(",")
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        filter_with_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${meta.id[tumor_idx]} > ${group}.agg.pon.vcf
        """
	stub:
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		"""
		echo ${meta.id[tumor_idx]}
		touch ${group}.agg.pon.vcf
		"""
        
}

process FFPE_PON_FILTER {
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"

	when:
		params.assay == "solid"

	input:
		tuple val(group), val(meta), file(vcf)
		
	output:
		tuple val(group), val(meta), file("${group}.agg.pon.ponffpe.vcf"), emit: vcf_pon_ffpe

	script:
        def pons = []
        if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
        if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
        if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
        def pons_str = pons.join(",")
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        filter_with_ffpe_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${meta.id[tumor_idx]} > ${group}.agg.pon.ponffpe.vcf
        """
	stub:
		tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
		"""
		echo ${meta.id[tumor_idx]}
		touch ${group}.agg.pon.ponffpe.vcf
		"""
        
}

process ANNOTATE_VEP {
	container = params.vepcon
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	label "process_high"
	tag "$group"
    
	input:
		tuple val(group), val(meta), file(vcf)
    
	output:
		tuple val(group), val(meta), file("${out}"), emit: vcf_vep

	script:
		out = vcf.getBaseName()
		out = out + ".vep.vcf"
		"""
		vep -i ${vcf} -o ${out} \\
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
		${params.custom_vep} \\
	"""
	stub:
		out = vcf.getBaseName()
		out = out + ".vep.vcf"
		"""
		echo ${params.custom_vep} $params.VEP_CACHE
		touch ${out}
		"""
}

process MARK_GERMLINES {
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	label "process_low"
	tag "$group"

	input:
		tuple val(group), val(meta), file(vcf) // from vcf_germline.join(meta_germline.groupTuple())

		
	output:
		tuple val(group), val(meta), file("${group}.agg.pon.vep.markgerm.vcf"), emit: vcf_germline


	script:
		if( meta.id.size() >= 2 ) {
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			fix_vep_gnomad.pl $vcf > ${group}.agg.pon.vep.fix.vcf
			mark_germlines.pl --vcf ${group}.agg.pon.vep.fix.vcf --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]} --assay $params.assay > ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			fix_vep_gnomad.pl $vcf > ${group}.agg.pon.vep.fix.vcf
			mark_germlines.pl --vcf ${group}.agg.pon.vep.fix.vcf --tumor-id ${meta.id[0]} --assay $params.assay > ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
	stub:
		if( meta.id.size() >= 2 ) {
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			echo --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]}
			touch ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
		else if( meta.id.size() == 1 ) {
			"""
			echo ${meta.id[0]}
			touch ${group}.agg.pon.vep.markgerm.vcf
			"""
		}
}

process FILTER_FOR_CNV {
    publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"

    input:
		tuple val(group), val(meta), file(vcf), val(vc), file(vcf_unfilt)

	output:
		tuple val(group), file("${group}_vardict.germlines.vcf.gz"), file("${group}_vardict.germlines.vcf.gz.tbi"), emit: vcf_only_germline
    
    script:
        """
        germline_for_cnvkit.pl $vcf > ${group}.agg.pon.vep.germline.vcf
		bedtools intersect -a $vcf_unfilt -b ${group}.agg.pon.vep.germline.vcf -header > ${group}_vardict.germlines.vcf
        bgzip ${group}_vardict.germlines.vcf
        tabix ${group}_vardict.germlines.vcf.gz
        """
	stub:
		"""
		echo $vcf $vcf_unfilt
		touch ${group}_vardict.germlines.vcf.gz ${group}_vardict.germlines.vcf.gz.tbi
		"""

}

process COYOTE_SEGMENTS {
	publishDir "${params.outdir}/${params.subdir}/cnv", mode: 'copy', overwrite: true
	label "process_single"
	tag "${meta.id}"

	input:
		tuple val(group), val(meta), file(vcf)
	
	output:
		tuple val(group), val(meta), file("${meta.id}.cn-segments.panel.bed"), emit: filtered
		tuple val(group), val(meta), file("${meta.id}.cn-segments.bed"), emit: raw

	script:
		normal = ""
		if ( meta.type == 'normal' || meta.type == 'N'  ) {
			normal = "--normal"
		}
		"""
		coyote_segmentator.pl --vcf $vcf --panel $params.panel_cna --id ${meta.id} $normal --genes /fs1/resources/ref/hg38/gtf/gencode.v33.annotation.genes.proteincoding.bed
		"""
	stub:
		panel = params.cnv_panel_path + "/" + meta.diagnosis + ".cna"
		"""
		touch ${meta.id}.cn-segments.panel.bed ${meta.id}.cn-segments.bed
		"""
}

process MERGE_SEGMENTS {
	publishDir "${params.outdir}/${params.subdir}/cnv", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(meta), file(segments)

	output:
		tuple val(group), file("${group}.cn-segments.panel.merged.bed"), emit: merged

	script:
		"""
		cat $segments > ${group}.cn-segments.panel.merged.bed
		"""
	stub:
		"""
		touch ${group}.cn-segments.panel.merged.bed
		"""

}

process FILTER_MANTA {
	publishDir "${params.outdir}/${params.subdir}/svvcf", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(meta), file(vcf)

	output:
		tuple val(group), val(meta), file("${meta.id}_manta_filtered.vcf"), emit: filtered
		tuple val(group), val(meta), file("${meta.id}_manta_bnd_filtered.vcf"), emit: bnd_filtered

	script:
		"""
		filter_manta.pl --vcf $vcf --id ${meta.id} --af 0.05
		"""

	stub:
		"""
		touch ${meta.id}_manta_bnd_filtered.vcf ${meta.id}_manta_filtered.vcf
		"""

}

process GENEFUSE_JSON_TO_VCF {
	publishDir "${params.outdir}/${params.subdir}/svvcf", mode: 'copy', overwrite: true
	label "process_single"
	tag "$group"
	
	input:
		tuple val(group), val(meta), file(json)

	output:
		tuple val(group), file("${meta.id}_genefuse.vcf"), emit:genefuse_vcf

	script:
		"""
		python /fs1/viktor/SomaticPanelPipeline_dsl2/bin/genefuse_json_to_vcf.py -i ${meta.id} -j $json -o ${meta.id}_genefuse.vcf
		"""

	stub:
		"""
		touch ${meta.id}_genefuse.vcf
		"""
}

process BIOMARKERS_TO_JSON {
	label "process_single"
	tag "$group"
	publishDir "${params.outdir}/${params.subdir}/biomarkers", mode: 'copy', overwrite: true

	input:
		tuple val(group), file(markers)

	output:
		tuple val(group), file("${group}.bio.json"), emit: biomarkers_json

	script:
		msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
		tmp = (msis ?: []) + (msip ?: []) + (hrd ?: [])
        command = tmp.join(' ')
		"""
		python /fs1/viktor/SomaticPanelPipeline_dsl2/bin/aggregate_biomarkers.py $command --out ${group}.bio.json --id $group
		"""
	stub:
		msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
		tmp = [msis, msip, hrd]
		tmp = tmp - null
        command = tmp.join(' ')
		"""
		echo $command
		touch ${group}.bio.json
		"""
}

process VCFANNO {
	label "process_single"
	tag "$group"


	input:
		tuple val(group), val(meta), file(vcf) 
		
	output:
		tuple val(group), val(meta), file("${group}.agg.enigma.vcf"), emit: vcf_enigma

	script:
		"""
		vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno $vcf > ${group}.agg.enigma.vcf
		"""
	stub:
		"""
		touch ${group}.agg.enigma.vcf
		"""
}

process CREATE_SNVPON {
	label "process_single"
	tag "$vc"


	input:
		tuple val(group), val(vc), file(vcfs) 
		
	output:
		tuple val(group), val(vc), file("${params.assay}_${vc}_PON.snv"), emit: SNV_PON

	script:
		"""
		create_snv_pon.pl "*.vcf.gz" > ${params.assay}_${vc}_PON.snv
		"""
	stub:
		"""
		echo $vcfs
		touch ${params.assay}_${vc}_PON.snv
		"""
}

process CONTAMINATION {
	publishDir "${params.outdir}/${params.subdir}/QC/contamination", mode: 'copy', overwrite: true, pattern: "*.png"
	publishDir "${params.outdir}/${params.subdir}/QC/contamination", mode: 'copy', overwrite: true, pattern: "*.txt"
	publishDir "${params.crondir}/contamination", mode: 'copy', overwrite: true, pattern: "*.contamination"
	container = "/fs1/resources/containers/perl-gd.sif"
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(meta), file(vcf)

	output:
		tuple val(group), file("*.txt"), file("*.png"), emit: contamination_result_files
		tuple val(group), file("*.contamination"), emit: contamination_cdm

	script:
		if(meta.id.size() >= 2) { 
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} --assay ${params.cdm} --detect-level 0.01 > ${meta.id[tumor_idx]}.value
			echo "--overwrite --sample-id ${meta.id[tumor_idx]} --run-folder ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --contamination" > ${meta.id[tumor_idx]}.1
			paste -d " " ${meta.id[tumor_idx]}.1 ${meta.id[tumor_idx]}.value > ${meta.id[tumor_idx]}.contamination
			find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} --assay ${params.cdm} --detect-level 0.01 --normal > ${meta.id[normal_idx]}.value
			echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --contamination" > ${meta.id[normal_idx]}.1
			paste -d " " ${meta.id[normal_idx]}.1 ${meta.id[normal_idx]}.value > ${meta.id[normal_idx]}.contamination
			"""
		}
		else {
			"""
			find_contaminant.pl --vcf $vcf --case-id ${meta.id[0]} --assay ${params.cdm} --detect-level 0.01 > ${meta.id[0]}.value
			echo "--overwrite --sample-id ${meta.id[0]} --sequencing-run ${meta.sequencing_run[0]} --assay ${params.cdm} --contamination" > ${meta.id[0]}.1
			paste -d " " ${meta.id[0]}.1 ${meta.id[0]}.value > ${meta.id[0]}.contamination
			"""
		}
	stub:
		if(meta.id.size() >= 2) { 
			tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
			"""
			touch test.png
			touch test.txt
			touch ${meta.id[tumor_idx]}.contamination
			touch ${meta.id[normal_idx]}.contamination
			"""
		}
		else {
			"""
			touch test.png
			touch test.txt
			touch ${meta.id[0]}.contamination
			"""
		}	


}