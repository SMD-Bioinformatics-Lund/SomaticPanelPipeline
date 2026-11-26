process PURECN_FORMAT_VCF {
    label "process_single"
    tag "${group}"
    
    input:
        tuple val(group), val(meta), path(vcf) // from process MARK_GERMLINES in SNV_ANNOTATE subworkflow
    
    output:
        tuple val(group), val(meta), path("*.agg.pon.vep.markgerm.purecn.vcf.gz"), path("*.agg.pon.vep.markgerm.purecn.vcf.gz.tbi"),    emit: vcf_purecn
        path "versions.yml",                                                                                                           emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${group}"
            if (meta.id.size() >= 2) {
                prefix = "${prefix}p"
                }
        // old name of the script: 5b_format_vcf_for_purecn_AD_POPAF.py.
        """
        format_vcf_for_purecn.py \\
            --input ${vcf} \\
            --output ${prefix}.agg.pon.vep.markgerm.ad_popaf.vcf
        
        bgzip -c ${prefix}.agg.pon.vep.markgerm.ad_popaf.vcf > ${prefix}.agg.pon.vep.markgerm.ad_popaf.vcf.gz
        
        bcftools index -t ${prefix}.agg.pon.vep.markgerm.ad_popaf.vcf.gz
        
        bcftools annotate \\
            --annotation ${params.annot_dbsnp} \\
            --columns ID \\
            --output ${prefix}.agg.pon.vep.markgerm.purecn.vcf.gz \\
            --output-type z \\
            ${prefix}.agg.pon.vep.markgerm.ad_popaf.vcf.gz
        
        bcftools index \\
            -t ${prefix}.agg.pon.vep.markgerm.purecn.vcf.gz
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
            bgzip: \$(bgzip --version | grep 'bgzip' | sed 's/.* //g')
            bcftools: \$(bcftools --version | head -n1 | sed -e 's/bcftools //')
        END_VERSIONS

        """


    stub:
        def prefix = task.ext.prefix ?: "${group}"
                if (meta.id.size() >= 2) {
                    prefix = "${prefix}p"
                    }
    
        """
        touch ${prefix}.agg.pon.vep.markgerm.purecn.vcf.gz
        touch ${prefix}.agg.pon.vep.markgerm.purecn.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
            bgzip: \$(bgzip --version | grep 'bgzip' | sed 's/.* //g')
            bcftools: \$(bcftools --version | head -n1 | sed -e 's/bcftools //')
        END_VERSIONS
        """
        
}


process PURECN_COVERAGE {
    label "process_single"
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), path(tumor_bam), path(tumor_bai)
    
    output:
        tuple val(group), val(meta), path(params.purecn_skip_gc_norm
            ? "${meta.id}.*_coverage.txt.gz"
            : "${meta.id}.*_coverage_loess.txt.gz"),                         emit: tumor_cov
        path "versions.yml",                                                 emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def gc_norm_option = params.purecn_skip_gc_norm ? '--skip-gc-norm' : ''
    
        """
        Rscript /usr/local/lib/R/site-library/PureCN/extdata/Coverage.R \\
            --bam ${tumor_bam} \\
            --intervals ${params.purecn_intervals} ${gc_norm_option}
        
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            purecn: \$(Rscript -e 'cat(as.character(packageVersion("PureCN")))')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.T.bwa.umi.sort_coverage${params.purecn_skip_gc_norm ? '' : '_loess'}.txt.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            purecn: \$(Rscript -e 'cat(as.character(packageVersion("PureCN")))')
        END_VERSIONS
        """

}


process PURECN_CALL {
    label "process_single"
    tag "${meta.id}"
    
    input:
        tuple val(group), path(purecn_vcf), path(purecn_tbi), val(meta), path(tumor_cov)
    
    output:
        tuple val(group), val(meta), path("${meta.id}_loh.csv"),                            emit: purecn_loh_raw_calls
        tuple val(group), val(meta), path("${meta.id}_purecn_snv_posteriors.tsv"),          emit: purecn_snv_posteriors
        tuple val(group), val(meta), path("${meta.id}.log"),                                emit: purecn_log
        tuple val(group), val(meta), path("${meta.id}.rds"),                                emit: purecn_rds
        path "versions.yml",                                                                emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''

        """
        Rscript /usr/local/lib/R/site-library/PureCN/extdata/PureCN.R \\
            --out . \\
            --tumor ${tumor_cov} \\
            --vcf ${purecn_vcf} \\
            --sampleid ${meta.id} \\
            --intervals ${params.purecn_intervals} \\
            --normaldb ${params.purecn_normaldb} \\
            --mapping-bias-file ${params.purecn_mapping_bias} \\
            $args
        
        Rscript -e '
        args <- commandArgs(trailingOnly=TRUE)
        res <- readRDS(args[1])
        id <- 1
        bm <- res\$results[[id]]\$SNV.posterior\$posteriors
        write.table(
            bm, 
            args[2],
            sep="\\t",
            quote=FALSE,
            row.names=FALSE
            )
        ' \\
        ${meta.id}.rds \\
        ${meta.id}_purecn_snv_posteriors.tsv

        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            purecn: \$(Rscript -e 'cat(as.character(packageVersion("PureCN")))')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_loh.csv
        touch ${meta.id}_purecn_snv_posteriors.tsv
        touch ${meta.id}.log
        touch ${meta.id}.rds
        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            purecn: \$(Rscript -e 'cat(as.character(packageVersion("PureCN")))')
        END_VERSIONS
        """

}
 
process PURECN_LOH_EXTRACT {
    label "process_single"
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), path(purecn_loh_raw_calls), path(purecn_log), path(purecn_rds)
    
    output:
        tuple val(group), val(meta), file("${meta.id}_purecn_loh_filtered.csv"),            emit: purecn_loh_filtered
        path "versions.yml",                                                                emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def snp_threshold = meta.id.size() >= 2 ? 10 : 10
        """
        filter_purecn_loh_output.py \\
            -l ${purecn_log} \\
            -p ${purecn_rds} \\
            -i ${purecn_loh_raw_calls} \\
            --snp_threshold ${snp_threshold} \\
            -o ${meta.id}_purecn_loh_filtered.csv

        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_purecn_loh_filtered.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

// needs to be run only if purecn_loh_filtered is not empty
process PURECN_CNLOH_GENS {
    label "process_single"
    tag "${meta.id}"
    
    input:
        tuple val(group), val(meta), path(purecn_loh_filtered), path(purecn_snv_posteriors)
    
    output:
        tuple val(group), val(meta), file("${meta.id}_informative_cnloh_snp.tsv"),                          emit: purecn_cnloh_snp
        tuple val(group), val(meta), file("${meta.id}_purecn_cnloh_gens.tsv"),                              emit: purecn_cnloh_gens
        tuple val(group), val(meta), file("${meta.id}_baf_purecn_cnlohseg.png"),                            emit: purecn_cnloh_baf_plot, optional: true
        tuple val(group), val(meta), file("${meta.id}_cnloh_coyote.json"),                                  emit: purecn_cnloh_coyote_json
        path "versions.yml",                                                                                emit: versions

    // checking the input contains data
    when:
        (task.ext.when == null || task.ext.when) && purecn_loh_filtered.text.readLines().size() > 1
    
    script: // old name of the script: 8c_CNLOH_extract_infoSNP_make_GENSbed_plotVAF.py. Change it also in the description of the python script
        """
        cnloh_visualization_and_extractSNP.py \\
            -i ${purecn_loh_filtered} \\
            -i2 ${purecn_snv_posteriors} \\
            --out_snp ${meta.id}_informative_cnloh_snp.tsv \\
            --out_gens ${meta.id}_purecn_cnloh_gens.tsv \\
            --out_coyote ${meta.id}_cnloh_coyote.json
 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS

        """

    stub:
        """
        touch ${meta.id}_informative_cnloh_snp.tsv
        touch ${meta.id}_purecn_cnloh_gens.tsv
        touch ${meta.id}_baf_purecn_cnlohseg.png
        touch ${meta.id}_cnloh_coyote.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed -e 's/Python //g')
        END_VERSIONS
        """
}

    
    
    


