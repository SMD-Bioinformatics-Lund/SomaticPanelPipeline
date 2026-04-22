
process PEDDY {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf), val(ped)

    output:
        tuple val(group), val(meta), path("${group}.ped_check.csv"),path("${group}.peddy.ped"), path("${group}.sex_check.csv"), emit: peddy_files

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        peddy --sites hg38 $vcf $ped --prefix ${prefix}
        """
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo $meta > ${prefix}.ped_check.csv
        echo $meta > ${prefix}.peddy.ped
        echo $meta > ${prefix}.sex_check.csv
        """
}