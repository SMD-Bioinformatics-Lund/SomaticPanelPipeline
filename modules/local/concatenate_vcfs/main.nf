process AGGREGATE_VCFS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcfs)

    output:
        tuple val(group), val(meta), file("*.agg.unsorted.vcf"),    emit: vcf_agg
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"

        sample_order = meta.id[0]
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
        }

        """
        aggregate_vcf.pl --vcfs ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} > ${prefix}.agg.unsorted.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            sample_order = meta.id[tumor_idx]+","+meta.id[normal_idx]
            """
            echo tumor:${meta.id[tumor_idx]} normal:${meta.id[normal_idx]}
            touch ${prefix}.agg.unsorted.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else {
            """
            touch ${prefix}.agg.unsorted.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}
