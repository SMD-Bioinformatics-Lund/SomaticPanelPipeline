process CALL_LOH {
    label "process_single"
    tag "${meta.id}"

    input: 
        // tumor-only: filtered to meta.type == "T" before calling this process
        tuple val(group), val(meta), val(part), file(cns)

    output:
        tuple val(group), val(meta), val(part), file("*.loh_cat.tsv"),  emit: loh_cat
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${group}.${meta.id}${part}"
        
            """
            cnvkit_segment_to_loh.py \
                --input ${cns} \
                --sample ${meta.id} \
                --sex ${meta.sex} \
                --out_cat ${prefix}.loh_cat.tsv

        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
    
    stub:
        def prefix = task.ext.prefix ?: "${group}.${meta.id}${part}"
        """
        touch ${prefix}.loh_cat.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        
        """
}