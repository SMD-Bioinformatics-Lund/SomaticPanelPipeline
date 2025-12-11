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
        tuple val(group), val(meta), file("*.${part}.baf.bed.gz"), file("*.${part}.cov.bed.gz"),    emit: cnvkit_gens
        path "versions.yml",                                                                        emit: versions
        
    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        generate_gens_data_from_cnvkit.pl $cnr $vcf ${meta.id}
        mv ${meta.id}.baf.bed.gz ${prefix}.${part}.baf.bed.gz
        mv ${meta.id}.cov.bed.gz ${prefix}.${part}.cov.bed.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.${part}.baf.bed.gz 
        touch ${prefix}.${part}.cov.bed.gz 
        touch ${prefix}.gens

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process CNVKIT_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(part), file(cns), file(cnr), file(vcf), file(tbi)
        val (tc)

    output:
        tuple val(group), val(meta), val(part), file("*.${part}.call*.cns"),            emit: cnvkitsegment
        tuple val(group), val(meta), val(part), file("*.${part}_logr_ballele.cnvkit"),  emit: cnvkit_baflogr
        tuple val(group), val(meta), val(part), file("*.${part}.cnvkit.vcf"),           emit: cnvkit_vcf
        path "versions.yml",                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"
        def prefix2  = task.ext.prefix2 ?: "${meta.id}.${meta.type}"

        call = "cnvkit.py call $cns -v $vcf -o ${prefix}.${part}.call.cns"
        callvcf = "cnvkit.py export vcf ${prefix}.${part}.call.cns -i '${meta.id}' > ${prefix2}.${part}.cnvkit.vcf"
        if (meta.purity && tc == "true") {
            call = "cnvkit.py call $cns -v $vcf --purity ${meta.purity} -o ${prefix}.${part}.call.purity.cns"
            callvcf = "cnvkit.py export vcf ${prefix}.${part}.call.purity.cns -i '${meta.id}' > ${prefix2}.${part}.cnvkit.vcf"
        }

        """
        $call
        $callvcf
        cnvkit.py export nexus-ogt -o ${prefix}.${part}_logr_ballele.cnvkit ${cnr} ${vcf}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix   = task.ext.prefix ?: "${group}.${meta.id}"
        def prefix2  = task.ext.prefix2 ?: "${meta.id}.${meta.type}"

        call = "cnvkit.py call $cns -v $vcf -o ${prefix}.${part}.call.cns \\ cnvkit.py export vcf ${prefix}.${part}.call.cns -i '${meta.id}' > ${prefix2}.${part}.vcf"
        if (meta.purity && tc == "true") {
            call = "cnvkit.py call $cns -v $vcf --purity ${meta.purity} -o ${prefix}.${part}.call.purity.cns \\ cnvkit.py export vcf ${prefix}.${part}.call.purity.cns -i '${meta.id}' > ${prefix2}.${part}.vcf"
        }

        """
        echo $call
        touch ${prefix}.${part}.call.purity.cns
        touch ${prefix}.${part}_logr_ballele.cnvkit 
        touch ${prefix2}.${part}.cnvkit.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process MERGE_GENS {
    label 'process_single'
    // THIS PROCESS SUCKS, please fix.
    input:
        tuple val(group), val(meta), file(baf), file(cov)

    output:
        tuple val(group), val(meta), file("*baf.bed.gz*"), file("*cov.bed.gz*"), optional: true,    emit: merged_gens
        tuple val(group), val(meta), file("*.gens"),                                                emit: dbload
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        process_group = group
        if ( meta.paired ) {
            process_group = group + 'p'
        }
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${meta.id}"

    shell:
        '''
        if !{params.cnvkit_split} == true
        then
        for i in $( ls *.cov.bed.gz ); do zgrep "^o_" $i | sed 's/o_//' >> !{meta.id}.base.cov.bed ; done
            bedtools sort -i !{meta.id}.base.cov.bed > !{meta.id}.base.cov.bed.sort
            sed 's/^/o_/' !{meta.id}.base.cov.bed.sort >> !{meta.id}.merged.cov.bed
            sed 's/^/a_/' !{meta.id}.base.cov.bed.sort >> !{meta.id}.merged.cov.bed
            sed 's/^/b_/' !{meta.id}.base.cov.bed.sort >> !{meta.id}.merged.cov.bed
            sed 's/^/c_/' !{meta.id}.base.cov.bed.sort >> !{meta.id}.merged.cov.bed
            sed 's/^/d_/' !{meta.id}.base.cov.bed.sort >> !{meta.id}.merged.cov.bed
            bedtools sort -i !{meta.id}.merged.cov.bed > !{meta.id}.merged.sorted.cov.bed
            bgzip !{meta.id}.merged.sorted.cov.bed
            tabix !{meta.id}.merged.sorted.cov.bed.gz
            for i in $( ls *.baf.bed.gz ); do zgrep "^o_" $i | sed 's/o_//' >> !{meta.id}.base.baf.bed ; done
            bedtools sort -i !{meta.id}.base.baf.bed > !{meta.id}.base.baf.bed.sort
            sed 's/^/o_/' !{meta.id}.base.baf.bed.sort >> !{meta.id}.merged.baf.bed
            sed 's/^/a_/' !{meta.id}.base.baf.bed.sort >> !{meta.id}.merged.baf.bed
            sed 's/^/b_/' !{meta.id}.base.baf.bed.sort >> !{meta.id}.merged.baf.bed
            sed 's/^/c_/' !{meta.id}.base.baf.bed.sort >> !{meta.id}.merged.baf.bed
            sed 's/^/d_/' !{meta.id}.base.baf.bed.sort >> !{meta.id}.merged.baf.bed
            bedtools sort -i !{meta.id}.merged.baf.bed > !{meta.id}.merged.sorted.baf.bed
            bgzip !{meta.id}.merged.sorted.baf.bed
            tabix !{meta.id}.merged.sorted.baf.bed.gz
            echo "gens load sample --sample-id !{meta.id} --case-id !{process_group} --genome-build 38 --baf !{params.gens_accessdir}/!{meta.id}.merged.sorted.baf.bed.gz --coverage !{params.gens_accessdir}/!{meta.id}.merged.sorted.cov.bed.gz" > !{meta.id}.gens
        else
            tabix !{meta.id}.full.baf.bed.gz
            tabix !{meta.id}.full.cov.bed.gz
            echo "gens load sample --sample-id !{meta.id} --case-id !{process_group} --genome-build 38 --baf !{params.gens_accessdir}/!{meta.id}.full.baf.bed.gz --coverage !{params.gens_accessdir}/!{meta.id}.full.cov.bed.gz" > !{meta.id}.gens
        fi

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            bedtools: $(bedtools --version | sed -e "s/bedtools v//g")
            bgzip: $(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: $(echo $(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*$//')
        END_VERSIONS
        '''

    stub:
        process_group = group
        if ( meta.paired ) {
            process_group = group + 'p'
        }
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${meta.id}"
        """
        echo "gens load sample --sample-id ${meta.id} --case-id ${process_group} --genome-build 38 --baf ${meta.id}.merged.sorted.baf.bed.gz --coverage ${meta.id}.merged.sorted.cov.bed.gz" > ${meta.id}.gens
        touch ${meta.id}.merged.sorted.cov.bed.gz 
        touch ${meta.id}.merged.sorted.baf.bed.gz 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}


process GENS_V4 {
    label 'process_single'
    // THIS PROCESS SUCKS, please fix.
    input:
        tuple val(group), val(meta), file(baf), file(cov)

    output:
        tuple val(group), val(meta), file("*.gens_v4"),     emit: dbload_v4

    when:
        task.ext.when == null || task.ext.when
    
    script:
        process_group = group
        if ( meta.paired ) {
            process_group = group + 'p'
        }
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${meta.id}"

        """
        echo "gens load sample --sample-id ${meta.id} --case-id ${process_group} --genome-build 38 --sample-type ${meta.type} --baf ${params.gens_accessdir}/${meta.id}.merged.sorted.baf.bed.gz --coverage ${params.gens_accessdir}/${meta.id}.merged.sorted.cov.bed.gz" > ${meta.id}.gens_v4
        """
    stub:
        process_group = group
        if ( meta.paired ) {
            process_group = group + 'p'
        }
        def args     = task.ext.args ?: ""
        def prefix   = task.ext.prefix ?: "${meta.id}"
        """
        echo "gens load sample --sample-id ${meta.id} --case-id ${process_group} --genome-build 38 --sample-type ${meta.type} --baf ${params.gens_accessdir}/${meta.id}.merged.sorted.baf.bed.gz --coverage ${params.gens_accessdir}/${meta.id}.merged.sorted.cov.bed.gz" > ${meta.id}.gens_v4
        """
}