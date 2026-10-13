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
        def caller_list = callers instanceof List ? callers : [callers]
        def vcf_list = vcfs instanceof List ? vcfs : [vcfs]
        def records = [caller_list, vcf_list].transpose()
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

process SVDB_QUERY {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), path(vcf)
        tuple val(meta2), path(db)
        
    output:
        tuple val(group), val(meta), path("out/*.vcf"), emit: queried_vcf
        path "versions.yml",                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'query.vcf'
        """
        mkdir -p out
        svdb --query ${args} \\
            --db ${db} \\
            --query_vcf ${vcf} \\
            > out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( svdb 2>&1 | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def suffix = task.ext.suffix ?: 'query.vcf'
        """
        mkdir -p out
        touch out/${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$( svdb 2>&1 | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        END_VERSIONS
        """
}
