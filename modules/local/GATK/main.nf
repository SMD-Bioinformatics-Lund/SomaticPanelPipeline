// SOMATIC CALLING //
process GATKCOV_BAF {
    cpus 10
    memory '64 GB'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("${meta.id}.allelicCounts.tsv"),  emit: gatk_baf
        path "versions.yml",                                                emit: versions
    
    script:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx50g" CollectAllelicCounts \\
            -L $params.GATK_GNOMAD \\
            -I $bam \\
            -R $params.genome_file \\
            -O ${meta.id}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${meta.id}.allelicCounts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATKCOV_COUNT {
    cpus 10
    memory '64 GB'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true 

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("${meta.id}.standardizedCR.tsv"), file("${meta.id}.denoisedCR.tsv"),  emit: gatk_count
        path "versions.yml",                                                                                    emit: versions

    script:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk CollectReadCounts \\
            -I ${bam} -L $params.gatk_intervals_full \\
            --interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5
        gatk --java-options '-Xmx50g' DenoiseReadCounts \\
            -I ${bam}.hdf5 --count-panel-of-normals $params.GATK_pon \\
            --standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
            --denoised-copy-ratios ${meta.id}.denoisedCR.tsv
        gatk PlotDenoisedCopyRatios \\
            --standardized-copy-ratios ${meta.id}.standardizedCR.tsv \\
            --denoised-copy-ratios ${meta.id}.denoisedCR.tsv \\
            --sequence-dictionary $params.GENOMEDICT \\
            --minimum-contig-length 46709983 --output . --output-prefix ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${meta.id}.standardizedCR.tsv ${meta.id}.denoisedCR.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """    
}


process GATKCOV_CALL {
    publishDir "${params.outdir}/${params.subdir}/gatkcov", mode: 'copy', overwrite: true, pattern: '*.seg'
    publishDir "${params.outdir}/${params.subdir}/plots", mode: 'copy', overwrite: true, pattern: '*.png'
    cpus 10
    memory '64 GB'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

    input:
        tuple val(group), val(meta), file(allele), file(stdCR), file(denoised)

    output:
        tuple val(group), file("${meta.id[tumor_idx]}.called.seg"),     emit: gatcov_called
        tuple val(group), file("${meta.id[tumor_idx]}.modeled.png"),    emit: gatcov_plot
        path "versions.yml",                                            emit: versions

    script:
        normalcounts = ""
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        if (meta.id.size() == 2) {
            normalcounts = "--normal-allelic-counts " + allele[normal_idx]
        }

        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx40g" ModelSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${allele[tumor_idx]} \\
            --minimum-total-allele-count-normal 20 \\
            --output . \\
            --output-prefix ${meta.id[tumor_idx]} \\
            $normalcounts
        gatk CallCopyRatioSegments \\
            --input ${meta.id[tumor_idx]}.cr.seg \\
            --output ${meta.id[tumor_idx]}.called.seg
        gatk PlotModeledSegments \\
            --denoised-copy-ratios ${denoised[tumor_idx]} \\
            --allelic-counts ${meta.id[tumor_idx]}.hets.tsv \\
            --segments ${meta.id[tumor_idx]}.modelFinal.seg \\
            --sequence-dictionary $params.GENOMEDICT \\
            --minimum-contig-length 46709983 \\
            --output . \\
            --output-prefix ${meta.id[tumor_idx]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T'  }
        normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N'  }
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${meta.id[tumor_idx]}.called.seg
        touch ${meta.id[tumor_idx]}.modeled.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


// GERMLINE CALLING //
process GATK_COUNT_GERMLINE {
    cpus 10
    memory '50GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("${meta.id}.tsv"),    emit: count_germline
        path "versions.yml",                                    emit: versions

    script:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx20g" CollectReadCounts \\
            -L $params.gatk_intervals \\
            -R $params.genome_file \\
            -imr OVERLAPPING_ONLY \\
            -I $bam \\
            --format TSV -O ${meta.id}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${meta.id}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process GATK_CALL_PLOIDY {
    cpus 10
    memory '40GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    tag "$group"

    input:
        tuple val(group), val(meta), file(tsv)

    output:
        tuple val(group), val(meta), file("ploidy.tar"),    emit: gatk_ploidy
        path "versions.yml",                                emit: versions

    script:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk
        gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
            --model $params.ploidymodel \\
            -I $tsv \\
            -O ploidy/ \\
            --output-prefix $group
        tar -cvf ploidy.tar ploidy/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ploidy.tar

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process GATK_CALL_GERMLINE_CNV {
    cpus 8
    memory '45GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: true, pattern: '*.vcf'
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(tsv), file(ploidy), val(i), val(refpart)

    output:
        tuple val(group), val(meta), val(i), file("${group}_${i}.tar"), emit: gatk_call_germline
        path "versions.yml",                                            emit: versions

    script:
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        tar -xvf ploidy.tar
        mkdir ${group}_${i}
        gatk --java-options "-Xmx25g" GermlineCNVCaller \\
            --run-mode CASE \\
            -I $tsv \\
            --contig-ploidy-calls ploidy/${group}-calls/ \\
            --model ${refpart} \\
            --output ${group}_${i}/ \\
            --output-prefix ${group}_${i}
        tar -cvf ${group}_${i}.tar ${group}_${i}/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        set +u
        source activate gatk
        export HOME=/local/scratch
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}

        touch ${group}_${i}.tar

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process POSTPROCESS {
    cpus 8
    memory '40GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: 'true'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), val(i), file(tar), file(ploidy), val(shard_no), val(shard)

    output:
        tuple val(group), val(meta), file("${meta.id}_gatk_genotyped-intervals.vcf.gz"),    emit: gatk_germline_intervals
        tuple val(group), val(meta), file("${meta.id}_gatk_genotyped-segments.vcf.gz"),     emit: gatk_germline_segmentsvcf
        tuple val(group), val(meta), file("${meta.id}_gatk_denoised.vcf.gz"),               emit: gatk_germline_denoised_log2
        path "versions.yml",                                                                emit: versions

    script:
        modelshards = shard.join(' --model-shard-path ') // join each reference shard
        caseshards = []
        for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
            tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls' 
            caseshards = caseshards + tmp
        }
        caseshards = caseshards.join( ' --calls-shard-path ')

    shell:
        '''
        THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
        for model in !{tar}; do
        tar -xvf $model
        done
        tar -xvf !{ploidy}
        set +u
        source activate gatk
        export MKL_NUM_THREADS=!{task.cpus}
        export OMP_NUM_THREADS=!{task.cpus}
        gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
            --allosomal-contig X --allosomal-contig Y \
            --contig-ploidy-calls ploidy/!{group}-calls/ \
            --sample-index 0 \\
            --output-genotyped-intervals !{meta.id}_gatk_genotyped-intervals.vcf.gz \
            --output-genotyped-segments !{meta.id}_gatk_genotyped-segments.vcf.gz \
            --output-denoised-copy-ratios !{meta.id}_gatk_denoised.vcf.gz \
            --sequence-dictionary !{params.GENOMEDICT} \
            --calls-shard-path !{caseshards} \
            --model-shard-path !{modelshards}

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            gatk4: $(echo $(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*$//')
            tar: $(echo $(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        '''

    stub:
        """
        export THEANO_FLAGS="base_compiledir=."
        export MKL_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        set +u
        source activate gatk

        touch ${meta.id}_gatk_genotyped-intervals.vcf.gz 
        touch ${meta.id}_gatk_genotyped-segments.vcf.gz 
        touch ${meta.id}_gatk_denoised.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            tar: \$(echo \$(tar --version) |sed 's/tar (GNU tar) //; s/ .*//')
        END_VERSIONS
        """
}


process FILTER_MERGE_GATK {
    publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: 'true'
    tag "${meta.id}"
    cpus 2
    memory '1GB'
    time '1h'

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), file("${meta.id}.gatk.filtered.merged.vcf"),  emit: gatk_normal_vcf
        path "versions.yml",                                            emit: versions

    script:
        """
        filter_gatk.pl $gatk > ${meta.id}.gatk.filtered.vcf
        mergeGATK.pl ${meta.id}.gatk.filtered.vcf > ${meta.id}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.gatk.filtered.merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

}


process GATK2VCF {
    publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: 'true'
    tag "${meta.id}"
    cpus 2
    memory '1GB'
    time '1h'

    input:
        tuple val(group), file(seg), val(meta)

    output:
        tuple val(group), val(meta), file("${meta.id}_gatk_tumor.vcf"), emit: tumor_vcf
        path "versions.yml",                                            emit: versions

    script:
        """
        python3 /fs1/viktor/SomaticPanelPipeline_dsl2/bin/gatk_to_vcf.py -s $seg -o ${meta.id}_gatk_tumor.vcf -p ${meta.purity}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_gatk_tumor.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process MERGE_GATK_TUMOR {
    publishDir "${params.outdir}/${params.subdir}/svvcf/", mode: 'copy', overwrite: 'true'
    tag "${meta.id}"
    cpus 2
    memory '1GB'
    time '1h'

    input:
        tuple val(group), val(meta), file(gatk)

    output:
        tuple val(group), val(meta), file("${meta.id}_gatk_tumor_merged.vcf"),  emit: tumor_vcf_merged
        path "versions.yml",                                                    emit: versions

    script:
        """
        mergeGATK_tumor.pl $gatk > ${meta.id}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_gatk_tumor_merged.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}