process SVDB_MERGE_PANEL {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcfs)

    output:
        tuple val(group), val(meta), file("*.merged.vcf"), emit: merged_vcf
        path "versions.yml",                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def prefix  = task.ext.prefix   ?: "${meta.id}"

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
        svdb --merge --vcf $vcfs_svdb $args --priority $priority > ${prefix}.merged.vcf
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
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
        echo $vcfs_svdb $priority
        touch ${prefix}.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

}


process SVDB_MERGE_SINGLES {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(vc), file(vcfs)
        
    output:
        tuple val(group), val(vc), file("*_${vc}.merged.vcf"),  emit: singles_merged_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}" 
        vcfs_svdb = vcfs.join(' ')
        vc = vc[0]
        """
        svdb --merge --vcf $vcfs_svdb $args > ${prefix}_${vc}.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        vcfs_svdb = vcfs.join(' ')
        vc = vc[0]
        """
        touch ${prefix}_${vc}.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

}

process SVDB_ANNOTATE_ARTEFACTS {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        
    output:
        tuple val(group), val(meta), file("${meta.id}.cnv.artefacts.vcf"),  emit: artefacts
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}" 
        """
        svdb --query $args --query_vcf $vcf > ${meta.id}.cnv.artefacts.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "-db $params.loqusdb_export --query_vcf $vcf" > ${meta.id}.cnv.artefacts.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

}