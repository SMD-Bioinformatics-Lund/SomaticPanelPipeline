process CNVKITREF {
    label 'process_high_memory'
    label 'process_medium_cpus'

    input:
        tuple val(name), val(id), file(bam), file(bai)
        val(part)
        path(bedfile)

    output:
        tuple val(name), file("*_cnvkit_${part}.cnn"),  emit: cnvkit_ref
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${name}"
        """
        cnvkit.py batch --normal *.bam \\
            --targets ${bedfile} \\
            $args \\
            --output-reference ${prefix}_cnvkit_${part}.cnn \\
            --output-dir results/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${name}"
        """
        echo $bedfile $part
        touch ${prefix}_cnvkit_${part}.cnn

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process CNVKIT_BATCH {
    label "process_low"
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)
        val(reference)
        val(part)

    output:
        tuple val(group), val(meta), file("*.${part}.cns"), val(part),  emit: cnvkit_cns
        tuple val(group), val(meta), file("*.${part}.cnr"), val(part),  emit: cnvkit_cnr
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${group}.${meta.id}"

        """
        cnvkit.py batch $bam -r $reference -d results/
        cp results/*sort.cnr ${prefix}.${part}.cnr
        cp results/*sort.cns ${prefix}.${part}.cns

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${group}.${meta.id}"
        """
        touch ${prefix}.${part}.cns
        touch ${prefix}.${part}.cnr

        echo $reference

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process CNVKIT_PLOT {
    label 'process_single'
    label "stage"
    label "scratch"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(cns), file(cnr), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.cnvkit_overview.png"),  emit: cnvkitplot
        path "versions.yml",                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${group}.${meta.id}"

        """
        cnvkit.py scatter -s *.cn{s,r} -o ${prefix}.${part}.cnvkit_overview.png -v ${vcf} -i ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${group}.${meta.id}"

        """
        touch ${prefix}.${part}.cnvkit_overview.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process CNVKIT_GENS {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(cnr), val(part), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.baf.bed"), file("*.${part}.cov.bed"), emit: cnvkit_gens
        path "versions.yml",                                                                         emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def out_prefix = part ? "${prefix}.${part}" : "${prefix}"
        """
        generate_gens_data_from_cnvkit_bed.pl --cnr $cnr --vcf $vcf --sample-id ${meta.id} --baf-out ${out_prefix}.baf.bed --cov-out ${out_prefix}.cov.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def out_prefix = part ? "${prefix}.${part}" : "${prefix}"
        """
        touch ${out_prefix}.baf.bed
        touch ${out_prefix}.cov.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process CNVKIT_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(cns), file(vcf), file(tbi)
        val (tc)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.call*.cns"),    emit: cnvkitsegment
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"
        call = "cnvkit.py call $cns -v $vcf -o ${prefix}.${part}.call.cns"
        if (tc == "true" && meta.purity) {
            call = "cnvkit.py call $cns -v $vcf --purity ${meta.purity} -o ${prefix}.${part}.call.purity.cns"
        }

        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        $call

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"
        def purity_suffix = tc == "true" && meta.purity ? '.purity' : ''

        """
        touch ${prefix}.${part}.call${purity_suffix}.cns

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process CNVKIT_EXPORT_VCF {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(cns)

    output:
        tuple val(group), val(meta), val(part), file("*.cnvkit*.vcf"),  emit: cnvkit_vcf
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args     ?: ""
        def prefix      = task.ext.prefix   ?: "${meta.id}.${meta.type}"
        def suffix      = task.ext.suffix   ?: "${part}.cnvkit.vcf"

        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        cnvkit.py export vcf ${cns} -i '${meta.id}' ${args} > ${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def args        = task.ext.args     ?: ""
        def prefix      = task.ext.prefix   ?: "${meta.id}.${meta.type}"
        def suffix      = task.ext.suffix   ?: "${part}.cnvkit.vcf"

        """
        echo "cnvkit.py export vcf ${cns} -i '${meta.id}' ${args} > ${prefix}.${suffix}"
        touch ${prefix}.${suffix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process CNVKIT_EXPORT_NEXUS_OGT {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(cnr), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}_logr_ballele.cnvkit"),  emit: cnvkit_baflogr, optional: true
        path "versions.yml",                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"

        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        if zcat -f ${vcf} | grep -qv '^#'; then
            cnvkit.py export nexus-ogt \
                -o ${prefix}.${part}_logr_ballele.cnvkit \
                ${cnr} \
                ${vcf}
        else
            echo "Skipping CNVkit nexus-ogt export: empty VCF"
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"

        """
        touch ${prefix}.${part}_logr_ballele.cnvkit

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}
