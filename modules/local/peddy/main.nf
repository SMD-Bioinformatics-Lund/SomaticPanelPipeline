
process PEDDY {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), path("*.ped_check.csv"),path("*.peddy.ped"), path("*.sex_check.csv"), emit: peddy_files
        path "versions.yml",                                                                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        def sex_code = (
            meta.sex == 'M' ? 1 :
            meta.sex == 'F' ? 2 : 0
        )
        def ped_sex = (
            meta.sex == 'M' ? 'male' :
            meta.sex == 'F' ? 'female' : 'unknown'
        )

        """
        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        echo -e "${group}\t${meta.id}\t0\t0\t${sex_code}\t2" > ${meta.id}.ped
        peddy --sites hg38 $vcf ${meta.id}.ped --prefix ${prefix}

        # Targeted panels can contain no Peddy sex-chromosome sites. Preserve the
        # completed Peddy reports and let the CDM conversion record sex QC as unavailable.
        if [[ ! -f ${prefix}.sex_check.csv ]]; then
            printf 'sample_id,ped_sex,predicted_sex,error\\n%s,%s,unknown,unavailable\\n' \\
                '${meta.id}' '${ped_sex}' > ${prefix}.sex_check.csv
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
        END_VERSIONS
        """
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo $meta > ${prefix}.ped_check.csv
        echo $meta > ${prefix}.peddy.ped
        echo $meta > ${prefix}.sex_check.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
        END_VERSIONS
        """
}
