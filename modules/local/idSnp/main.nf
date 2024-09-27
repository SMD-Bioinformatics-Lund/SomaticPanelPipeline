process ALLELE_CALL {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*final.vcf"),    emit:   sample_id_vcf
        path "versions.yml",                                emit:   versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""
        def args2   = task.ext.args2 ?: ""
        def args3   = task.ext.args3 ?: ""
        """
        bcftools mpileup $args $bam | bcftools call $args2   > ${prefix}.raw.vcf
        bcftools annotate $args3 -o ${prefix}.final.vcf ${prefix}.raw.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.final.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/bcftools //; s/ .*//')
        END_VERSIONS
        """
}

process SNP_CHECK {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcfs)

    output:
        tuple val(group), val(meta), file("*.json"),            emit: idsnp_checked
        path "versions.yml", optional: true,                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        
        def args    = task.ext.args  ?: ""
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
	    normalvcf   = vcfs[normal_idx]
	    tumorvcf    = vcfs[tumor_idx]

        if(meta.id.size() == 2) {
            """
            idsnp_controller-myeolid_specific.pl \\
                --vcf_sample $tumorvcf  \\
                --vcf_control $normalvcf \\
                --sample  $tumor_id \\
                --control $normal_id \\
                $args

            cp s${tumor_id}_c${normal_id}.json  ${tumor_id}.json
            cp s${tumor_id}_c${normal_id}.json  ${normal_id}.json
            rm s${tumor_id}_c${normal_id}.json
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  "{ "Info" :  "False"}" > ${tumor_id}.json
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
        
        if(meta.id.size() == 2) {
            """ 
            touch s${tumor_id}_c${normal_id}.csv
            touch  ${tumor_id}.json
            touch  ${normal_id}.json


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        } else {
        
            """
            echo "Not applicable" > s${tumor_id}.csv
            echo  "{ "Info" :  "False"}" > ${tumor_id}.json

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}


process PAIR_CDM {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(jsons)

    output:
        tuple val(group), val(meta), file("*.pairscore"), emit: isnsp_cdm_done

    when:
        task.ext.when == null || task.ext.when

    script:
      def args    = task.ext.args  ?: ""
      if(meta.id.size() == 2) {
        
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
        normal_id   = meta.id[normal_idx]
        tumor_id    = meta.id[tumor_idx]
	    normaljson  = jsons[normal_idx]
	    tumorjson    = jsons[tumor_idx]

        """
        echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${tumorjson}" > ${meta.id[tumor_idx]}.pairscore
            
        echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${normaljson}" > ${meta.id[normal_idx]}.pairscore
        """

      } else {
        tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        tumor_id    = meta.id[tumor_idx]
	    tumorjson    = jsons[tumor_idx]
        """
        echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${tumorjson}" > ${meta.id[tumor_idx]}.pairscore
        """
      }
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        
        if(meta.id.size() == 2) {
            tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx  = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            normal_id   = meta.id[normal_idx]
            tumor_id    = meta.id[tumor_idx]
	        normaljson  = jsons[normal_idx]
	        tumorjson    = jsons[tumor_idx]
             """
             echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${tumorjson}" > ${meta.id[tumor_idx]}.pairscore
             
             echo "--overwrite --sample-id ${meta.id[normal_idx]} --sequencing-run ${meta.sequencing_run[normal_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${normaljson}" > ${meta.id[normal_idx]}.pairscore
            """
        } else {
            tumor_idx   = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            tumor_id    = meta.id[tumor_idx]
	        tumorjson    = jsons[tumor_idx]
            """
            echo "--overwrite --sample-id ${meta.id[tumor_idx]} --sequencing-run ${meta.sequencing_run[tumor_idx]} --assay ${params.cdm} --pairedness ${params.outdir}/${params.subdir}/QC/${tumorjson}" > ${meta.id[tumor_idx]}.pairscore
            """
        }

}

process PROVIDER {
    label 'process_single'
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    
    output:
        tuple val(group), val(meta), file("*.genotypes"),    emit: genotype_checked
        path "versions.yml",                                 emit: versions  // Emit version information in YAML format
    
    when:
        task.ext.when == null || task.ext.when
        
    script:
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args  ?: ""

        """
        provider.pl  --bam $bam  $args  --out $prefix 
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
        
    stub:
        def prefix  = task.ext.prefix ?: "${meta.id}"

        """
        touch ${prefix}.genotypes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}