process COYOTE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(importy)

    output:
        tuple val(group), file("*.coyote"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        normal_sample = null
        sample_no = meta.id.size()
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
        }
        tumor_sample = meta.id[tumor_idx]

        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio' }
        lowcov_idx     = importy.findIndexOf{ it =~ 'lowcov' }
        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvkit_idx     = importy.findIndexOf{ it =~ 'cnvkit_overview.png' }

        // add matching flags //
        cnvseg     = cnvseg_idx     >= 0 ? importy[cnvseg_idx].collect {'--cnv /access/' + params.subdir + '/cnv/' + it } : null
        fusions    = fusions_idx    >= 0 ? importy[fusions_idx].collect {'--transloc /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'--biomarkers /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot    = cnvplot_idx    >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null
        lowcov     = lowcov_idx     >= 0 ? importy[lowcov_idx].collect {'--lowcov /access/' + params.subdir + '/QC/' + it } : null
        purity     = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { '--purity ' + it} : null
        tmp        = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null

        if ( cnvplot == null && cnvkit_idx >= 0 ){
            cnvplot = importy[cnvkit_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it }
        }

        lowcov = lowcov_idx >= 0 ? importy[lowcov_idx].collect {'--lowcov /access/' + params.subdir + '/QC/' + it } : null
        purity = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { '--purity ' + it} : null
        tmp = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        import_command = tmp.join(' ')

        """
        echo "/data/bnf/scripts/import_DSL2_to_coyote.pl --group $params.coyote_group \\
            --vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
            --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
            --build 38 \\
            --gens ${meta.id[tumor_idx]} \\
            --subpanel ${meta.diagnosis[tumor_idx]} \\
            --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} \\
            --sample_no $sample_no \\
            --case_id ${tumor_sample} \\
            --control_id ${normal_sample} \\
            --profile ${environment} \\
            --assay $params.coyote_group \\
            $import_command" > ${process_group}.coyote
        """

    stub:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        normal_sample = null
        sample_no = meta.id.size()
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
        }
        tumor_sample = meta.id[tumor_idx]
        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx     = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx    = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio' }
        lowcov_idx     = importy.findIndexOf{ it =~ 'lowcov' }
        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        cnvkit_idx     = importy.findIndexOf{ it =~ 'cnvkit_overview.png' }


        // add matching flags //
        cnvseg     = cnvseg_idx     >= 0 ? importy[cnvseg_idx].collect {'--cnv /access/' + params.subdir + '/cnv/' + it } : null
        fusions    = fusions_idx    >= 0 ? importy[fusions_idx].collect {'--transloc /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'--biomarkers /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot    = cnvplot_idx    >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null
        lowcov     = lowcov_idx     >= 0 ? importy[lowcov_idx].collect {'--lowcov /access/' + params.subdir + '/QC/' + it } : null
        purity     = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { '--purity ' + it} : null
        tmp        = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null

        if ( cnvplot == null && cnvkit_idx >= 0 ){
            cnvplot = importy[cnvkit_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it }
        } 

        lowcov = lowcov_idx >= 0 ? importy[lowcov_idx].collect {'--lowcov /access/' + params.subdir + '/QC/' + it } : null
        purity = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { '--purity ' + it} : null
        tmp = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        import_command = tmp.join(' ')

        """        
        echo "/data/bnf/scripts/import_DSL2_to_coyote.pl --group $params.coyote_group \\
            --vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
            --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
            --build 38 \\
            --gens ${meta.id[tumor_idx]} \\
            --subpanel ${meta.diagnosis[tumor_idx]} \\
            --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} \\
            --sample_no $sample_no \\
            --case_id ${tumor_sample} \\
            --control_id ${normal_sample} \\
            --profile ${environment} \\
            --assay $params.coyote_group \\
            $import_command" > ${process_group}.coyote
        """
}

process COYOTE_YAML {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(importy)

    output:
        tuple val(group), file("*.coyote3.yaml"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        normal_sample = null
        sample_no = meta.id.size()
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
        }
        tumor_sample = meta.id[tumor_idx]
        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx     = importy.findIndexOf{ it =~ 'cnvs' }
        fusions_idx    = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio.json' }
        cov_idx        = importy.findIndexOf{ it =~ 'cov.json' }
        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        lowcov_idx     = importy.findIndexOf{ it =~ 'lowcov.bed' }
        cnvkit_idx     = importy.findIndexOf{ it =~ 'cnvkit_overview.png' }

        // add matching flags //
        cnvseg     = cnvseg_idx     >= 0 ? importy[cnvseg_idx].collect {'cnv: /access/' + params.subdir + '/cnv/' + it } : null
        fusions    = fusions_idx    >= 0 ? importy[fusions_idx].collect {'transloc: /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'biomarkers: /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot    = cnvplot_idx    >= 0 ? importy[cnvplot_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it } : null
        lowcov     = lowcov_idx     >= 0 ? importy[lowcov_idx].collect {'lowcov: /access/' + params.subdir + '/QC/' + it } : null
        cov        = cov_idx        >= 0 ? importy[cov_idx].collect {'cov: /access/' + params.subdir + '/QC/' + it } : null
        purity     = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { 'purity: ' + it} : null
        tmp        = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it } : null
        if ( cnvplot == null && cnvkit_idx >= 0 ){
            cnvplot = importy[cnvkit_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it }
        } 
        lowcov = lowcov_idx >= 0 ? importy[lowcov_idx].collect {'lowcov: /access/' + params.subdir + '/QC/' + it } : null
        purity = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { 'purity: ' + it} : null
        tmp = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [] ) + (cov ?: [] )
        import_command = tmp.join('\n')

        """
        echo --- > ${process_group}.coyote3.yaml
        echo groups: [\\'$params.coyote_group\\'] >> ${process_group}.coyote3.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote3.yaml
        echo clarity-sample-id: \\'${meta.clarity_sample_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo clarity-pool-id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo genome_build: 38 >> ${process_group}.coyote3.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote3.yaml
        echo sample_no: ${sample_no} >> ${process_group}.coyote3.yaml
        echo case_id: \\'${tumor_sample}\\' >> ${process_group}.coyote3.yaml
        echo control_id: \\'${normal_sample}\\' >> ${process_group}.coyote3.yaml
        echo profile: \\'${environment}\\' >> ${process_group}.coyote3.yaml
        echo assay: \\'$params.coyote_group\\' >> ${process_group}.coyote3.yaml
        printf "$import_command" >> ${process_group}.coyote3.yaml
        """
    stub:
        environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        normal_sample = null
        sample_no = meta.id.size()
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_sample = meta.id[normal_idx]
            sample_no = meta.id.size()
        }
        tumor_sample = meta.id[tumor_idx]
        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx     = importy.findIndexOf{ it =~ 'cnvs' }
        fusions_idx    = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio.json' }
        cnvplot_idx    = importy.findIndexOf{ it =~ 'modeled.png' }
        lowcov_idx     = importy.findIndexOf{ it =~ 'lowcov.bed' }
        cov_idx        = importy.findIndexOf{ it =~ 'cov.json' }
        cnvkit_idx     = importy.findIndexOf{ it =~ 'cnvkit_overview.png' }

        // add matching flags //
        cnvseg     = cnvseg_idx     >= 0 ? importy[cnvseg_idx].collect {'cnv: /access/' + params.subdir + '/cnv/' + it } : null
        fusions    = fusions_idx    >= 0 ? importy[fusions_idx].collect {'transloc: /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'biomarkers: /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot    = cnvplot_idx    >= 0 ? importy[cnvplot_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it } : null
        lowcov     = lowcov_idx     >= 0 ? importy[lowcov_idx].collect {'lowcov: /access/' + params.subdir + '/QC/' + it } : null
        cov        = cov_idx        >= 0 ? importy[cov_idx].collect {'cov: /access/' + params.subdir + '/QC/' + it } : null
        purity     = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { 'purity: ' + it} : null
        tmp        = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it } : null
        if ( cnvplot == null && cnvkit_idx >= 0 ){
            cnvplot = importy[cnvkit_idx].collect {'cnvprofile: /access/' + params.subdir + '/plots/' + it }
        } 
        lowcov = lowcov_idx >= 0 ? importy[lowcov_idx].collect {'lowcov: /access/' + params.subdir + '/QC/' + it } : null
        purity = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { 'purity: ' + it} : null
        tmp = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [] ) + (cov ?: [] )
        import_command = tmp.join('\n')

        """
        echo --- > ${process_group}.coyote3.yaml
        echo groups: [\\'$params.coyote_group\\'] >> ${process_group}.coyote3.yaml
        echo subpanel: \\'${meta.diagnosis[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo name: \\'${process_group}\\' >> ${process_group}.coyote3.yaml
        echo clarity-sample-id: \\'${meta.clarity_sample_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo clarity-pool-id: \\'${meta.clarity_pool_id[tumor_idx]}\\' >> ${process_group}.coyote3.yaml
        echo genome_build: 38 >> ${process_group}.coyote3.yaml
        echo vcf_files: /access/${params.subdir}/vcf/${vcf} >> ${process_group}.coyote3.yaml
        echo sample_no: ${sample_no} >> ${process_group}.coyote3.yaml
        echo case_id: \\'${tumor_sample}\\' >> ${process_group}.coyote3.yaml
        echo control_id: \\'${normal_sample}\\' >> ${process_group}.coyote3.yaml
        echo profile: \\'${environment}\\' >> ${process_group}.coyote3.yaml
        echo assay: \\'$params.coyote_group\\' >> ${process_group}.coyote3.yaml
        printf "$import_command" >> ${process_group}.coyote3.yaml
        """
}
