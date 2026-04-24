
process PEDDY {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), path("*.ped_check.csv"),path("*.peddy.ped"), path("*.sex_check.csv"), emit: peddy_files

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        
        def sex_code = (
            meta.sex == 'M' ? 1 :
            meta.sex == 'F' ? 2 : 0
        )

        """
        echo -e "${group}\t${meta.id}\t0\t0\t${sex_code}\t2" > ${meta.id}.ped
        peddy --sites hg38 $vcf ${meta.id}.ped --prefix ${prefix}
        """
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo $meta > ${prefix}.ped_check.csv
        echo $meta > ${prefix}.peddy.ped
        echo $meta > ${prefix}.sex_check.csv
        """
}