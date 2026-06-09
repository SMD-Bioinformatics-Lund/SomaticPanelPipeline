process CALL_LOH {
    label "process_single"
    tag "$group"

    input: 
        tuple val(group), val(meta), val(part), file(cns)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.loh.bed"),  emit: loh_bed
        tuple val(group), val(meta), val(part), file("*.${part}.loh.txt"),  emit: loh_txt
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        
        """
        cnvkit_segments_loh.py \
            --input ${cns} \
            --sample ${meta.id} \
            --out ${group}.${meta.id}.${part}.loh.txt
            --out_bed ${group}.${meta.id}.${part}.loh.bed

        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
    
    stub:
        """
        touch ${group}.${meta.id}.${part}.loh.txt
        touch ${group}.${meta.id}.${part}.loh.bed
        
        """
}