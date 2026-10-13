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
        tuple val(group), val(meta), val(tumor_id), val(tumor_seq_run), file("${tumor_id}.T.idsnp"),    emit: idsnp_tumor, optional: true
        tuple val(group), val(meta), val(normal_id), val(normal_seq_run), file("${normal_id}.N.idsnp"), emit: idsnp_normal, optional: true
        path "versions.yml", optional: true,                                                            emit: versions

    when:
        (task.ext.when == null || task.ext.when) && SampleGroups.validateMetadata(group, meta)

    script:
        
        def args                    = task.ext.args  ?: ""
        def samples                 = SampleGroups.validateMetadata(group, meta)
        tumor_idx                   = samples.findIndexOf { SampleGroups.isTumor(it) }
        normal_idx                  = samples.findIndexOf { SampleGroups.isNormal(it) }
        normal_id                   = normal_idx >= 0 ? samples[normal_idx].id : null
        tumor_id                    = tumor_idx >= 0 ? samples[tumor_idx].id : null
                    
        tumor_seq_run               = tumor_idx >= 0 ? samples[tumor_idx].sequencing_run : null
        normal_seq_run              = normal_idx >= 0 ? samples[normal_idx].sequencing_run : null

        normalvcf                   = normal_idx >= 0 ? vcfs[normal_idx] : null
        tumorvcf                    = tumor_idx >= 0 ? vcfs[tumor_idx] : null

        norGTjson                   = normal_idx >= 0 ? genotypes[normal_idx] : null
        tumGTjson                   = tumor_idx >= 0 ? genotypes[tumor_idx] : null
        def single_id               = samples[0].id
        def single_role             = SampleGroups.isNormal(samples[0]) ? 'N' : 'T'

        if(samples.size() == 2) {
            """
            idsnp_controller.pl \\
                --vcf_sample $tumorvcf  \\
                --vcf_control $normalvcf \\
                --sample  $tumor_id \\
                --control $normal_id \\
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
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        } else {
            """
            echo "Not applicable" > s${single_id}.csv
            echo  '{ "is_paired_sample" : false}' > ${single_id}.json
            combinejsons.py --info-json ${single_id}.json --genotype-json ${genotypes[0]} --out ${single_id}.${single_role}.idsnp
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        def samples                 = SampleGroups.validateMetadata(group, meta)
        tumor_idx                   = samples.findIndexOf { SampleGroups.isTumor(it) }
        normal_idx                  = samples.findIndexOf { SampleGroups.isNormal(it) }
        normal_id                   = normal_idx >= 0 ? samples[normal_idx].id : null
        tumor_id                    = tumor_idx >= 0 ? samples[tumor_idx].id : null
        tumor_seq_run               = tumor_idx >= 0 ? samples[tumor_idx].sequencing_run : null
        normal_seq_run              = normal_idx >= 0 ? samples[normal_idx].sequencing_run : null
        def single_id               = samples[0].id
        def single_role             = SampleGroups.isNormal(samples[0]) ? 'N' : 'T'
        
        if(samples.size() == 2) {
            """ 
            touch s${tumor_id}_c${normal_id}.csv
            touch ${tumor_id}.T.idsnp
            touch ${normal_id}.N.idsnp


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        } else {
        
            """
            echo "Not applicable" > s${single_id}.csv
            echo  '{ "is_paired_sample" : false }' > ${single_id}.json
            touch ${single_id}.${single_role}.idsnp

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
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
