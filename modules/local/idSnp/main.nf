process GENOTYPE2JSON {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf), file(genotypes)

    output:
        tuple val(group), val(meta), file(vcf), file("*.genotypes.json"), emit: sample_id_genotypes
        path "versions.yml",                                             emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        """
        genotype2json.py --input ${genotypes} --out ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.genotypes.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process SNP_CHECK {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcfs), file(genotypes)

    output:
        tuple val(group), val(meta), val(tumor_id), val(tumor_seq_run), file("${tumor_id}.T.idsnp"),    emit: idsnp_tumor
        tuple val(group), val(meta), val(normal_id), val(normal_seq_run), file("${normal_id}.N.idsnp"), emit: idsnp_normal, optional: true
        path "versions.yml", optional: true,                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        
        def args                    = task.ext.args  ?: ""
        tumor_idx                   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx                  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id                   = meta.id[normal_idx]
        tumor_id                    = meta.id[tumor_idx]
                    
        tumor_seq_run               = meta.sequencing_run[tumor_idx]
        normal_seq_run              = meta.sequencing_run[normal_idx]

	    normalvcf                   = vcfs[normal_idx]
	    tumorvcf                    = vcfs[tumor_idx]

        norGTjson                   = genotypes[normal_idx]
        tumGTjson                   = genotypes[tumor_idx]

        if(meta.id.size() == 2) {
            """
            idsnp_controller.py \\
                --vcf-sample $tumorvcf  \\
                --vcf-control $normalvcf \\
                --sample  $tumor_id \\
                --control $normal_id \\
                --csv-out s${tumor_id}_c${normal_id}.csv \\
                --json-out s${tumor_id}_c${normal_id}.json \\
                $args

            ## Why is the perl-script not doing all of the below code?
            cp s${tumor_id}_c${normal_id}.json  ${tumor_id}.json
            cp s${tumor_id}_c${normal_id}.json  ${normal_id}.json
            rm s${tumor_id}_c${normal_id}.json
            today_date=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")

            echo '{"partner" : "${tumor_id}","sequencing_run" : "${tumor_seq_run}","analysis_date" : "'\${today_date}'" }' > ${normal_id}_partner_info.json
            echo '{"partner" : "${normal_id}","sequencing_run" : "${normal_seq_run}","analysis_date" : "'\${today_date}'" }' > ${tumor_id}_partner_info.json
            
            combinejsons.py --info-json ${tumor_id}.json --genotype-json ${tumGTjson} --partner-run-json-file  ${tumor_id}_partner_info.json --out ${tumor_id}.T.idsnp
            combinejsons.py --info-json ${normal_id}.json --genotype-json ${norGTjson} --partner-run-json-file ${normal_id}_partner_info.json --out ${normal_id}.N.idsnp
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        } else {
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  '{ "is_paired_sample" : false}' > ${tumor_id}.json
            combinejsons.py --info-json ${tumor_id}.json --genotype-json ${tumGTjson} --out ${tumor_id}.T.idsnp
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx                   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx                  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id                   = meta.id[normal_idx]
        tumor_id                    = meta.id[tumor_idx]
        tumor_seq_run               = meta.sequencing_run[tumor_idx]
        normal_seq_run              = meta.sequencing_run[normal_idx]
        
        if(meta.id.size() == 2) {
            """ 
            touch s${tumor_id}_c${normal_id}.csv
            touch ${tumor_id}.T.idsnp
            touch ${normal_id}.N.idsnp


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        } else {
        
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  '{ "is_paired_sample" : false }' > ${tumor_id}.json
            touch ${tumor_id}.T.idsnp

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
}

process PAIRGEN_CDM {
    label 'process_single'
    tag "${id}"

    input:
        tuple val(group), val(meta), val(id), val(run), file(json)

    output:
        tuple val(group), val(meta), file("${id}.pairgen"), emit: isnsp_cdm_done

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        echo "--overwrite --sample-id ${id} --sequencing-run ${run} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${json} " > ${id}.pairgen
        """

    stub:
        """
        echo "--overwrite --sample-id ${id} --sequencing-run ${run} --assay ${params.cdm} --id-snp ${params.outdir}/${params.subdir}/QC/${json} " > ${id}.pairgen
        """
}
