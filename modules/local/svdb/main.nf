process SVDB_MERGE_PANEL {
	label "process_single"
	tag "$group"
	publishDir "${params.outdir}/${params.subdir}/svvcf/merged/", mode: 'copy', overwrite: 'true'

	input:
		tuple val(group), val(meta), file(vcfs)

	output:
		tuple val(group), val(meta), file("${meta.id}.merged.vcf"), emit: merged_vcf

	script:
        def args = task.ext.args ?: ''  
        // for each sv-caller add idx, find vcf and find priority, add in priority order! //
        // index of vcfs added from mix //
        manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
        delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
        cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
        gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }
        genefuse_idx = vcfs.findIndexOf{ it =~ 'genefuse' }

        // find vcfs //
        manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
        delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
        cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
        gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
        genefuse = genefuse_idx >= 0 ? vcfs[genefuse_idx].collect {it + ':genefuse ' } : null
        tmp = (manta ?: []) + (delly ?: []) + (gatk ?: []) + (cnvkit ?: []) + (genefuse ?: [])
        vcfs_svdb = tmp.join(' ')

        // find priorities //
        mantap = manta_idx >= 0 ? 'manta' : null
        dellyp = delly_idx >= 0 ? 'delly' : null
        gatkp = gatk_idx >= 0 ? 'gatk' : null
        cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
        genefusep = genefuse_idx >= 0 ? 'genefuse' : null
        tmpp = [mantap, dellyp, gatkp, cnvkitp, genefusep]
        tmpp = tmpp - null
        priority = tmpp.join(',')
    
        """
        svdb --merge --vcf $vcfs_svdb --no_intra --pass_only $args --overlap 0.7 --priority $priority > ${meta.id}.merged.vcf
        """
    stub:
        manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
        delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
        cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
        gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }
        genefuse_idx = vcfs.findIndexOf{ it =~ 'genefuse' }

        // find vcfs //
        manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
        delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
        cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
        gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
        genefuse = genefuse_idx >= 0 ? vcfs[genefuse_idx].collect {it + ':genefuse ' } : null
        tmp = manta + delly + gatk + cnvkit
        tmp = tmp - null
        vcfs_svdb = tmp.join(' ')

        // find priorities //
        mantap = manta_idx >= 0 ? 'manta' : null
        dellyp = delly_idx >= 0 ? 'delly' : null
        gatkp = gatk_idx >= 0 ? 'gatk' : null
        cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
        genefusep = genefuse_idx >= 0 ? 'genefuse' : null
        tmpp = [mantap, dellyp, gatkp, cnvkitp, genefusep]
        tmpp = tmpp - null
        priority = tmpp.join(',')
        """
        echo $vcfs_svdb $priority
        touch ${meta.id}.merged.vcf
        """

}


process SVDB_MERGE_SINGLES {
	label "process_single"
	tag "$group"

	input:
		tuple val(group), val(vc), file(vcfs)
        
	output:
		tuple val(group), val(vc), file("${group}_${vc}.merged.vcf"), emit: singles_merged_vcf

	script:  
        vcfs_svdb = vcfs.join(' ')
        vc = vc[0]
        """
        svdb --merge --vcf $vcfs_svdb --no_intra --pass_only --bnd_distance 10 --overlap 1.0 > ${group}_${vc}.merged.vcf
        """
    stub:
        vcfs_svdb = vcfs.join(' ')
        vc = vc[0]
        """
        touch ${group}_${vc}.merged.vcf
        """

}