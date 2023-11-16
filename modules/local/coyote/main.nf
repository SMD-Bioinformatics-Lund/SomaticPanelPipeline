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
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        }
        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio' }
        cnvplot_idx = importy.findIndexOf{ it =~ 'png' }
        lowcov_idx = importy.findIndexOf{ it =~ 'lowcov' }

        // add matching flags //
        cnvseg = cnvseg_idx >= 0 ? importy[cnvseg_idx].collect {'--cnv /access/' + params.subdir + '/cnv/' + it } : null
        fusions = fusions_idx >= 0 ? importy[fusions_idx].collect {'--transloc /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'--biomarkers /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null
        lowcov = lowcov_idx >= 0 ? importy[lowcov_idx].collect {'--lowcov /access/' + params.subdir + '/QC/' + it } : null
        purity = meta.purity[tumor_idx] != false ? meta.purity[tumor_idx].toFloat().collect { '--purity ' + it} : null
        tmp = (cnvseg ?: []) + (fusions ?: []) + (biomarkers ?: []) + (cnvplot ?: []) + (lowcov ?: []) + (purity ?: [])
        import_command = tmp.join(' ')

        //echo "import_myeloid_to_coyote_vep_gms.pl --group $params.coyote_group \\
        """
        echo "/data/bnf/scripts/import_DSL2_to_coyote.pl --group $params.coyote_group \\
            --vcf /access/${params.subdir}/vcf/${vcf} --id ${process_group} \\
            --clarity-sample-id ${meta.clarity_sample_id[tumor_idx]} \\
            --build 38 \\
            --gens ${meta.id[tumor_idx]} \\
            --subpanel ${meta.diagnosis[tumor_idx]} \\
            --clarity-pool-id ${meta.clarity_pool_id[tumor_idx]} \\
            $import_command" > ${process_group}.coyote
        """

    stub:
        process_group = group
        tumor_idx = 0
        tumor_idx_lowcov = 0
        if( meta.id.size() >= 2 ) {
            process_group = group + 'p'
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            }
        // find what to load into coyote, depending on what files are in $import //
        // index of imports added from mix //
        cnvseg_idx = importy.findIndexOf{ it =~ 'panel' }
        fusions_idx = importy.findIndexOf{ it =~ 'annotated' }
        biomarkers_idx = importy.findIndexOf{ it =~ 'bio' }
        cnvplot_idx = importy.findIndexOf{ it =~ 'png' }
        lowcov_idx = importy.findIndexOf{ it =~ 'lowcov' }

        // add matching flags //
        cnvseg = cnvseg_idx >= 0 ? importy[cnvseg_idx].collect {'--cnv /access/' + params.subdir + '/cnv/' + it } : null
        fusions = fusions_idx >= 0 ? importy[fusions_idx].collect {'--transloc /access/' + params.subdir + '/fusions/' + it } : null
        biomarkers = biomarkers_idx >= 0 ? importy[biomarkers_idx].collect {'--biomarkers /access/' + params.subdir + '/biomarkers/' + it } : null
        cnvplot = cnvplot_idx >= 0 ? importy[cnvplot_idx].collect {'--cnvprofile  /access/' + params.subdir + '/plots/' + it } : null
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
            $import_command" > ${process_group}.coyote
        """
}