process CNVKIT2OVAHRDSCAR {
    label "process_single"
    tag "$id"

    input:
        tuple val(group), val(id), val(type), file(segments)

    output:
        tuple val(group), val(id), val(type), val(caller), file("*.cnvkit.ovaHRDscar.txt"), emit: ovaHRDscar_segments
        path "versions.yml",                                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${id}"

        caller = "cnvkit"
        if (segments =~ /purity/) {
            caller = "cnvkitpurity"
        }

        """
        cnvkit2HRD.pl $segments $id > ${prefix}.cnvkit.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
    
    stub:
        def prefix = task.ext.prefix ?: "${id}"
        caller = "cnvkit"
        if (segments =~ /purity/) {
            caller = "cnvkitpurity"
        }

        """
        touch ${prefix}.cnvkit.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}


process CNVKIT2SCARHRD {
    label "process_single"
    tag "$id"

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
        cnvkit2HRD.pl $segments ${meta.id} $ploidyv > ${prefix}.cnvkit.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
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
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
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

// TODO SHOULD BE REMOVED?
process OVAHRDSCAR {
    publishDir "${params.outdir}/${params.subdir}/ovaHRDscar/", mode: 'copy', overwrite: true, pattern: '*.txt'
    cpus 1
    time '1h'
    tag "$id"
    container = '/fs1/resources/containers/ovaHRDscar.sif'

    input:
        tuple val(group), val(id), val(type), val(sc), file(segments)

    output:
        tuple val(group), val(id), val(type), file("${group}_${sc}_ovaHRDscar_results.txt"),    emit: ovaHRDscar_score
        path "versions.yml",                                                                    emit: versions

    script:
        """
        Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/Run_ovarHRDscar.R $segments > ${group}_${sc}_ovaHRDscar_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}_${sc}_ovaHRDscar_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """
}