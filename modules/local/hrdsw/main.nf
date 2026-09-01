process CNVKIT2SCARHRD {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(segments)

    output:
        tuple val(group), val(meta), val(caller), file("${meta.id}.cnvkit.scarHRD.txt"),    emit: scarHRD_segments
        path "versions.yml",                                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        ploidyv = "NA"
        caller = "cnvkit"

        """
        cnvkit2HRD.py --input $segments --sample-id ${meta.id} --ploidy $ploidyv --out ${prefix}.cnvkit.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """  

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        ploidyv = "NA"
        caller = "cnvkit"
        """
        touch ${prefix}.cnvkit.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """ 
}


process SCARHRD {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(sc), file(segments)

    output:
        tuple val(group), file("*_scarHRD_results.txt"),    emit: scarHRD_score
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"

        """
        Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/Run_scarHRD.R $segments
        mv ${meta.id}_HRDresults.txt ${prefix}_${sc}_scarHRD_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_${sc}_scarHRD_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """ 

}
