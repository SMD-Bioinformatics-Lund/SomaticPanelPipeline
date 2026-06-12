process CALL_LOH {
    label "process_single"
    tag "$group"

    input: 
        tuple val(group), val(meta), val(part), file(cns)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.loh_cat.tsv"),  emit: loh_cat
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        
        """
        cnvkit_segments_loh.py \
            --input ${cns} \
            --sample ${meta.id} \
            --sex ${meta.sex} \
            --out_cat ${group}.${meta.id}.${part}.loh_cat.tsv

        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
    
    stub:
        """
        touch ${group}.${meta.id}.${part}.loh_cat.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        
        """
}