process VERIFYBAMID {
    label 'process_low'
    label 'scratch'
    label 'stage'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("${meta.id}.contaminationpy"),                emit: contamination
        tuple val(group), file("result.selfSM"), file("result.Ancestry"),               emit: results
        path "versions.yml",                                                            emit: versions
    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""   // reference 
        def args2   = task.ext.args2    ?: ""   // loci to check
        def args3   = task.ext.args3    ?: ""   // sanity check
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        verifybamid2 \
            $args3 \
            $args2 \
            $args \
            --BamFile $bam
        
        cut -f 7 result.selfSM | grep -v FREEMIX > contamination.value
        echo "--overwrite --sample-id ${meta.id} --sequencing-run ${meta.sequencing_run} --assay ${params.cdm} --contamination "> contamination.import
        paste contamination.import contamination.value > ${prefix}.contaminationpy
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid: \$(echo \$(verifybamid2 -h 2>&1 | grep Version | sed "s/Version://"))
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}.contaminationpy
        touch result.selfSM result.Ancestry
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid: \$(echo \$(verifybamid2 -h 2>&1 | grep Version | sed "s/Version://"))
        END_VERSIONS
        """
}