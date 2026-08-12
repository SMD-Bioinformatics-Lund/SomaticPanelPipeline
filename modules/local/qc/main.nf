
process QC_TO_CDM {
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

process SENTIEON_QC_TO_CDM {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(qc_files), path(cov_sample_summary), path(dedup)

    output:
        tuple val(group), val(meta), path("*_${meta.type}.QC"), emit: qc_cdm

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        qc_sentieon.py --sample-id ${meta.id}_${meta.type} --type panel --out ${prefix}_${meta.type}.QC
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_${meta.type}.QC
        """
}

process PANEL_DEPTH {
    label 'process_low'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(depth_tsv)

    output:
        tuple val(group), val(meta), file("*.lowcov.depth.bed"),    emit: lowcov_bed
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ""
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        panel_depth.py --depth $depth_tsv --out ${prefix}.lowcov.depth.bed $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.lowcov.depth.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """
}

process LOWCOV {
    label 'process_medium'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bed)

    output:
        tuple val(group), val(meta), val(meta.type), file("*.lowcov.bed"),  emit: lowcov_regions
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        overlapping_genes.py \\
            --intersect-bed $bed \\
            --out ${prefix}.lowcov.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.lowcov.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
        END_VERSIONS
        """
}

process QC_VALUES {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(qc)

    output:
        tuple val(group), val(meta), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV), emit: qc_melt_val
    
    when:
        task.ext.when == null || task.ext.when

    script:
        // Collect qc-data if possible from normal sample, if only tumor; tumor
        def ins_dev
        def coverage
        def ins_size
        qc.readLines().each{
            if (it =~ /\"(ins_size_dev)\"\s*:\s*\"(\S+)\"/) {
                ins_dev = it =~ /\"(ins_size_dev)\"\s*:\s*\"(\S+)\"/
            }
            if (it =~ /\"(mean_coverage)\"\s*:\s*\"(\S+)\"/) {
                coverage = it =~ /\"(mean_coverage)\"\s*:\s*\"(\S+)\"/
            }
            if (it =~ /\"(ins_size)\"\s*:\s*\"(\S+)\"/) {
                ins_size = it =~ /\"(ins_size)\"\s*:\s*\"(\S+)\"/
            }
        }
        INS_SIZE = ins_size[0][2]
        MEAN_DEPTH = coverage[0][2]
        COV_DEV = ins_dev[0][2]
        """
        echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
        """

    stub:
        // Collect qc-data if possible from normal sample, if only tumor; tumor
        def ins_dev
        def coverage
        def ins_size
        INS_SIZE = 1
        MEAN_DEPTH = 1
        COV_DEV = 1
        """
        echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
        """
}

process CONTAMINATION {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.txt"), file("*.png"),  emit: contamination_result_files
        tuple val(group), val(meta), file("*.value"),               emit: contamination_values
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''

        if(meta.id.size() >= 2) { 
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} $args > ${meta.id[tumor_idx]}.value
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[tumor_idx]} $args2 > ${meta.id[normal_idx]}.value

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            find_contaminant.pl --vcf $vcf --case-id ${meta.id[0]} $args > ${meta.id[0]}.value

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        if(meta.id.size() >= 2) { 
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            touch test.png
            touch test.txt
            touch ${meta.id[tumor_idx]}.value
            touch ${meta.id[normal_idx]}.value

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            touch test.png
            touch test.txt
            touch ${meta.id[0]}.value

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}

process CONTAMINATION_CDM {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(values)

    output:
        tuple val(group), val(meta), file("*.contaminationpy"), emit: contamination_cdm

    when:
        task.ext.when == null || task.ext.when

    script:
        if(meta.id.size() >= 2) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --contamination" > ${meta.id[tumor_idx]}.1
            paste -d " " ${meta.id[tumor_idx]}.1 ${meta.id[tumor_idx]}.value > ${meta.id[tumor_idx]}.contaminationpy

            echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --contamination" > ${meta.id[normal_idx]}.1
            paste -d " " ${meta.id[normal_idx]}.1 ${meta.id[normal_idx]}.value > ${meta.id[normal_idx]}.contaminationpy
            """
        }
        else {
            """
            echo "--overwrite --sample-id ${meta.id[0]} --sequencing-run ${meta.sequencing_run[0]} --assay ${params.cdm} --contamination" > ${meta.id[0]}.1
            paste -d " " ${meta.id[0]}.1 ${meta.id[0]}.value > ${meta.id[0]}.contaminationpy
            """
        }

    stub:
        if(meta.id.size() >= 2) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            touch ${meta.id[tumor_idx]}.contaminationpy
            touch ${meta.id[normal_idx]}.contaminationpy
            """
        }
        else {
            """
            touch ${meta.id[0]}.contaminationpy
            """
        }
}

process PEDDY2CDM {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), path(ped_check), path(peddy_ped), path(sex_check)

    output:
        tuple val(group), val(meta), path("*peddy.json"), emit: json
        tuple val(group), val(meta), path("*peddy2cdm"), emit: cdm

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''
        def args2   = task.ext.args2    ?: ''
        def sample_arg = "${meta.id}:${meta.sequencing_run}"
        """
        peddy2cdm.py \\
            --ped $ped_check \\
            --sex $sex_check \\
            --sample $sample_arg \\
            $args \\
            $args2
        """

    stub:
        """
        touch ${group}_peddy.json
        touch ${group}.peddy2cdm
        """
}

process LOWCOV_D4 {
    label 'process_medium'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), val(meta.type), file("*.cov.json"),    emit: coyote_cov_json
        tuple val(group), val(meta), val(meta.type), file("*.d4"),          emit: d4_coverage
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ""
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        coyote_d4_cov.py -b $bam $args -o ${prefix}.cov.json -s ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
            d4tools: \$(echo \$( d4tools 2>&1 | head -1 ) | sed "s/.*version: //" | sed "s/)//" )
            bedtools: \$(bedtools | grep Version | sed -r "s/Version:\s+//")
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix   ?: "${meta.id}"
        """
        touch ${prefix}.cov.json
        touch ${prefix}.d4

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
            d4tools: \$(echo \$( d4tools 2>&1 | head -1 ) | sed "s/.*version: //" | sed "s/)//" )
            bedtools: \$(bedtools | grep Version | sed -r "s/Version:\s+//")
        END_VERSIONS
        """
}
