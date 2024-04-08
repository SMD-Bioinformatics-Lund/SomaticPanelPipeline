process snpid {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(qc)

    output:
        tuple val(group), val(meta), file("*.cdmpy"), emit: cdm_done

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${prefix}.cdmpy
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "--sequencing-run ${meta.sequencing_run} --sample-type ${meta.type} --sample-id ${meta.id} --assay $params.cdm --qc ${params.outdir}/${params.subdir}/QC/$qc --lims-id ${meta.clarity_sample_id}" > ${prefix}.cdmpy
        """
}
