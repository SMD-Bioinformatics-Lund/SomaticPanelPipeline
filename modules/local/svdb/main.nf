process SVDB_MERGE_PANEL {
	cpus 1
	tag "$group"
	publishDir "${params.outdir}/${params.subdir}/svvcf/merged/", mode: 'copy', overwrite: 'true'
	time '10m'
	memory '1 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	input:
		tuple val(group), val(meta), file(vcfs)

	output:
		tuple val(group), val(meta), file("${meta.id}.merged.vcf")

	script:  
        // for each sv-caller add idx, find vcf and find priority, add in priority order! //
        // index of vcfs added from mix //
        manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
        delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
        cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
        gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

        // find vcfs //
        manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
        delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
        cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
        gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
        tmp = manta + delly + gatk + cnvkit
        tmp = tmp - null
        vcfs_svdb = tmp.join(' ')

        // find priorities //
        mantap = manta_idx >= 0 ? 'manta' : null
        dellyp = delly_idx >= 0 ? 'delly' : null
        gatkp = gatk_idx >= 0 ? 'gatk' : null
        cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
        tmpp = [mantap, dellyp, gatkp, cnvkitp]
        tmpp = tmpp - null
        priority = tmpp.join(',')
        // vcfs_svdb = cool
        // priority = cool
    
        """
        svdb --merge --vcf $vcfs_svdb --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority $priority > ${meta.id}.merged.vcf
        """
    stub:
        manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
        delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
        cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
        gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

        // find vcfs //
        manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
        delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
        cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
        gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
        tmp = manta + delly + gatk + cnvkit
        tmp = tmp - null
        vcfs_svdb = tmp.join(' ')

        // find priorities //
        mantap = manta_idx >= 0 ? 'manta' : null
        dellyp = delly_idx >= 0 ? 'delly' : null
        gatkp = gatk_idx >= 0 ? 'gatk' : null
        cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
        tmpp = [mantap, dellyp, gatkp, cnvkitp]
        tmpp = tmpp - null
        priority = tmpp.join(',')
        """
        echo $vcfs_svdb $priority
        touch ${meta.id}.merged.vcf
        """

}
