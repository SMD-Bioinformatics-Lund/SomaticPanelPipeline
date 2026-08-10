process SVDB_MERGE_PANEL {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), val(callers), file(vcfs)

    output:
        tuple val(group), val(meta), file("*.merged.vcf"), emit: merged_vcf
        path "versions.yml",                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        def caller_order = ['manta', 'delly', 'gatk', 'cnvkit', 'genefuse']
        def records = [callers, vcfs].transpose()
        def ordered_records = caller_order.collectMany { caller ->
            records.findAll { record -> record[0] == caller }
        }

        def vcfs_svdb = ordered_records.collect { record -> "${record[1]}:${record[0]}" }.join(' ')
        def priority = ordered_records.collect { record -> record[0] }.unique().join(',')

        """
        svdb --merge --vcf $vcfs_svdb $args --priority $priority > ${prefix}.merged.vcf
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
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
        tuple val(group), val(meta), val(vc), file(vcfs)
        
    output:
        tuple val(group), val(meta), val(vc), file("*_${vc}.merged.vcf"),   emit: singles_merged_vcf
        path "versions.yml",                                                emit: versions

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
        val(dbs)
        
    output:
        tuple val(group), val(meta), file("*.cnv.artefacts.vcf"),  emit: artefacts
        path "versions.yml",                                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def db_manifest = dbs.collect { db_name, db_vcf -> "${db_name}\\t${db_vcf}\\n" }.join('')

        """
        printf '${db_manifest}' > svdb_dbs.tsv

        svdb_annotate_artefacts.sh \\
            --query-vcf ${vcf} \\
            --dbs svdb_dbs.tsv \\
            --out ${prefix}.cnv.artefacts.vcf \\
            -- ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( svdb 2>&1 | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.cnv.artefacts.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( svdb 2>&1 | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """
}
