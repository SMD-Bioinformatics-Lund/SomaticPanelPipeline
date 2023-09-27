process CNVKIT2OVAHRDSCAR {
    cpus 1
    time '1h'
    tag "$id"

    input:
        tuple val(group), val(id), val(type), file(segments)

    output:
        tuple val(group), val(id), val(type), val(caller), file("${id}.cnvkit.ovaHRDscar.txt"), emit: ovaHRDscar_segments
        path "versions.yml",                                                                    emit: versions

    script:
        caller = "cnvkit"
        if (segments =~ /purity/) {
            caller = "cnvkitpurity"
        }

        """
        cnvkit2HRD.pl $segments $id > ${id}.cnvkit.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
    
    stub:
        caller = "cnvkit"
        if (segments =~ /purity/) {
            caller = "cnvkitpurity"
        }

        """
        touch ${id}.cnvkit.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}

process CNVKIT2SCARHRD {
    cpus 1
    time '1h'
    tag "$id"

    input:
        tuple val(group), val(meta), val(part), file(segments)

    output:
        tuple val(group), val(meta), val(caller), file("${meta.id}.cnvkit.scarHRD.txt"),    emit: scarHRD_segments
        path "versions.yml",                                                                emit: versions

    script:
        ploidyv = "NA"
        caller = "cnvkit"

        """
        cnvkit2HRD.pl $segments ${meta.id} $ploidyv > ${meta.id}.cnvkit.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """  

    stub:
        ploidyv = "NA"
        caller = "cnvkit"
        """
        touch ${meta.id}.cnvkit.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """ 
}

process ASCAT2SCARHRD {
    publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/segments", mode: 'copy', overwrite: true, pattern: '*.txt'
    cpus 1
    time '1h'
    tag "$id"

    input:
        tuple val(group), val(id), val(type), file(logr), file(baf) //val(ploidy)

    output:
        tuple val(group), val(id), val(type), val("ascat"), file("${id}.ascat.scarHRD.txt"),    emit: scarHRD_segments
        path "versions.yml",                                                                    emit: versions

    script:
        //ploidyv = ploidy.getText().trim()
        ploidyv = "NA"

        """
        ascat2HRD.pl $logr $baf $id $ploidyv > ${id}.ascat.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        ploidyv = "NA"

        """
        touch ${id}.ascat.scarHRD.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}

process ASCAT2OVAHRDSCAR {
    publishDir "${params.outdir}/${params.subdir}/ASCAT3.0/segments", mode: 'copy', overwrite: true, pattern: '*.txt'
    cpus 1
    time '1h'
    tag "$id"

    input:
        tuple val(group), val(id), val(type), file(logr), file(baf)

    output:
        tuple val(group), val(id), val(type), val("ascat"), file("${id}.ascat.ovaHRDscar.txt"), emit: ovaHRDscar_segments
        path "versions.yml",                                                                    emit: versions

    script:
        """
        ascat2HRD.pl $logr $baf $id > ${id}.ascat.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        """
        touch ${id}.ascat.ovaHRDscar.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}


process SCARHRD {
    publishDir "${params.outdir}/${params.subdir}/scarHRD/", mode: 'copy', overwrite: true, pattern: '*.txt'
    cpus 1
    time '1h'
    tag "${meta.id}"
    container = '/fs1/resources/containers/scarHRD.sif'

    input:
        tuple val(group), val(meta), val(sc), file(segments)

    output:
        tuple val(group), file("${meta.id}_${sc}_scarHRD_results.txt"), emit: scarHRD_score
        path "versions.yml",                                            emit: versions

    script:
        """
        Rscript /fs1/viktor/SomaticPanelPipeline_dsl2/bin/Run_scarHRD.R $segments
        mv ${meta.id}_HRDresults.txt ${meta.id}_${sc}_scarHRD_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_${sc}_scarHRD_results.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            Rscript: \$( Rscript --version |& sed 's/^.*version //' | cut -f1 -d ' ')
        END_VERSIONS
        """ 

}

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