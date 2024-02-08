process PHARMACOGENOMICS {
    label "process_single"
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bams), file(bais), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.pgx*.csv"), file("*.pgx*.sh"), emit: pgx_files

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args     ?: ''               // Anticipated bam output directory
        def args2   = task.ext.args2    ?: ''               // profiles to be used for analysis
        def args3   = task.ext.args3    ?: ''               // Anticipated csv output directory
        def suffix  = task.ext.suffix   ?: ".panel"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            tumor_bam = bams[tumor_idx].baseName
            normal_bam = bams[normal_idx].baseName
            prefix = "${meta.group[normal_idx]}"
            """
            echo -e "clarity_sample_id,id,type,assay,group,bam,bai,purity" > ${prefix}.pgx${suffix}.csv
            echo -e "${meta.clarity_sample_id[normal_idx]},${meta.id[normal_idx]},N,GMS_PGx_Panel,${meta.group[normal_idx]},${args}/${normal_bam}.bam,${args}/${normal_bam}.bam.bai,${meta.purity[normal_idx]}" >> ${prefix}.pgx${suffix}.csv

            echo -e "/fs2/sw/bnf-scripts/start_nextflow_analysis.pl $args3/${prefix}.pgx${suffix}.csv " > ${prefix}.pgx${suffix}.sh

            """
        }
        else if( meta.id.size() == 1 ) {
            bam     = bams[0].baseName
            prefix  = "${meta.group[0]}"
            """
            echo -e "clarity_sample_id,id,type,assay,group,bam,bai,purity" > ${prefix}.pgx${suffix}.csv
            echo -e "${meta.clarity_sample_id[0]},${meta.id[0]},${meta.type[0]},GMS_PGx_Panel,${meta.group[0]},${args}/${bam}.bam,${args}/${bam}.bam.bai,${meta.purity[0]}" >> ${prefix}.pgx${suffix}.csv

            echo -e "/fs2/sw/bnf-scripts/start_nextflow_analysis.pl $args3/${prefix}.pgx${suffix}.csv " > ${prefix}.pgx${suffix}.sh
            """
        }

    stub:
        def args    = task.ext.args     ?: ''               // Anticipated bam output directory
        def args2   = task.ext.args2    ?: ''               // profiles to be used for analysis
        def args3   = task.ext.args3    ?: ''               // Anticipated csv output directory
        def suffix  = task.ext.suffix   ?: ".panel"
        if( meta.id.size() >= 2 ) {
            tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            tumor_bam   = bams[tumor_idx].baseName
            normal_bam  = bams[normal_idx].baseName
            prefix      = "${meta.group[normal_idx]}"
            """
            echo -e "clarity_sample_id,id,type,assay,group,bam,bai,purity" > ${prefix}.pgx${suffix}.csv
            echo -e "${meta.clarity_sample_id[normal_idx]},${meta.id[normal_idx]},N,GMS_PGx_Panel,${meta.group[normal_idx]},${args}/${normal_bam}.bam,${args}/${normal_bam}.bam.bai,${meta.purity[normal_idx]}" >> ${prefix}.pgx${suffix}.csv

            echo -e "#/fs2/sw/bnf-scripts/start_nextflow_analysis.pl $args3/${prefix}.pgx${suffix}.csv " > ${prefix}.pgx${suffix}.sh
            """
        }
        else {
            bam     = bams[0].baseName
            prefix  = "${meta.group[0]}"
            """
            echo -e "clarity_sample_id,id,type,assay,group,bam,bai,purity" > ${prefix}.pgx${suffix}.csv
            echo -e "${meta.clarity_sample_id[0]},${meta.id[0]},${meta.type[0]},GMS_PGx_Panel,${meta.group[0]},${args}/${bam}.bam,${args}/${bam}.bam.bai,${meta.purity[0]}" >> ${prefix}.pgx${suffix}.csv

            echo -e "#/fs2/sw/bnf-scripts/start_nextflow_analysis.pl $args3/${prefix}.pgx${suffix}.csv " > ${prefix}.pgx${suffix}.sh
            """
        }
}