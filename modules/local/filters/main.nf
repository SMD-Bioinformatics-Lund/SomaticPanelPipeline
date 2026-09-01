process PON_FILTER {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) 
        
    output:
        tuple val(group), val(meta), file("*.agg.pon.vcf"), emit: vcf_pon
        path "versions.yml",                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        def pons = []
        if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
        if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
        if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
        def pons_str = pons.join(",")
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        filter_with_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${meta.id[tumor_idx]} > ${prefix}.agg.pon.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        echo ${meta.id[tumor_idx]}
        touch ${prefix}.agg.pon.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}


process FFPE_PON_FILTER {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        
    output:
        tuple val(group), val(meta), file("*.agg.pon.ponffpe.vcf"), emit: vcf_pon_ffpe
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        def pons = []
        if( params.freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
        if( params.vardict )   { pons.push("vardict="+params.PON_vardict) }
        if( params.tnscope )   { pons.push("tnscope="+params.PON_tnscope) }
        def pons_str = pons.join(",")
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        filter_with_ffpe_pon.pl --vcf $vcf --pons $pons_str --tumor-id ${meta.id[tumor_idx]} > ${prefix}.agg.pon.ponffpe.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
        """
        echo ${meta.id[tumor_idx]}
        touch ${prefix}.agg.pon.ponffpe.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """ 
}


process FIX_VEP_GNOMAD {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.fix.vcf"), emit: fixed_vcf
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        """
        fix_vep_gnomad.pl $vcf > ${prefix}.agg.pon.vep.fix.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.agg.pon.vep.fix.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process MARK_GERMLINES {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf) // from vcf_germline.join(meta_germline.groupTuple())

        
    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.markgerm.vcf"),    emit: vcf_germline
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            mark_germlines.pl --vcf $vcf --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]} $args > ${prefix}p.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            mark_germlines.pl --vcf $vcf --tumor-id ${meta.id[0]} $args > ${prefix}.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            echo --tumor-id ${meta.id[tumor_idx]} --normal-id ${meta.id[normal_idx]}
            touch ${prefix}p.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            echo ${meta.id[0]}
            touch ${prefix}.agg.pon.vep.markgerm.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
}


process GERMLINE_FOR_CNVKIT {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.agg.pon.vep.germline.vcf"), emit: germline_vcf
        path "versions.yml",                                             emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        """
        germline_for_cnvkit.pl $vcf > ${prefix}.agg.pon.vep.germline.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.agg.pon.vep.germline.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}

process FILTER_FREEBAYES {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), file("filtered_*.vcf"), emit: filtered_vcf
        path "versions.yml",                                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def part = vcf.name.replaceFirst(/\.filt1\.vcf$/, '')

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            filter_freebayes_somatic.py --vcf $vcf --tumor ${meta.id[tumor_idx]} --normal ${meta.id[normal_idx]} --out filtered_${part}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            filter_freebayes_unpaired.py --vcf $vcf --out filtered_${part}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        def part = vcf.name.replaceFirst(/\.filt1\.vcf$/, '')
        """
        touch filtered_${part}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process FILTER_VARDICT {
    label "process_low"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), file("vardict_*.vcf"), emit: vcfparts_vardict
        path "versions.yml",                                         emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        output_vcf = vcf.name.replaceFirst(/\.raw$/, '')

        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }

            """
            filter_vardict_somatic.pl $vcf ${meta.id[tumor_idx]} ${meta.id[normal_idx]} > ${output_vcf}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            filter_vardict_unpaired.pl $vcf > ${output_vcf}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                perl: \$(echo \$(perl -v 2>&1) | sed 's/.*(v//; s/).*//')
            END_VERSIONS
            """
        }

    stub:
        output_vcf = vcf.name.replaceFirst(/\.raw$/, '')
        """
        touch ${output_vcf}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process FILTER_PINDEL {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcf)

    output:
        tuple val(group), val(meta), val(vc), file("*_pindel.vcf"), emit: pindel_vcf
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        """
        filter_pindel_somatic.py --vcf ${vcf} --out ${prefix}_pindel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/.*\\s//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}_pindel.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/.*\\s//')
        END_VERSIONS
        """
}

process FILTER_TNSCOPE {
    label "process_medium"
    tag "$group"

    input:
        tuple val(group), val(meta), val("tnscope"), file(vcf)

    output:
        tuple val(group), val(meta), val("tnscope"), file("*.vcf"), emit: vcfparts_tnscope_filtered
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        if( meta.id.size() >= 2 ) {
            tumor_idx = meta.type.findIndexOf{ it == 'tumor' || it == 'T' }
            normal_idx = meta.type.findIndexOf{ it == 'normal' || it == 'N' }
            """
            filter_tnscope_somatic.py --vcf $vcf --tumor ${meta.id[tumor_idx]} --normal ${meta.id[normal_idx]} --out ${vcf}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }
        else if( meta.id.size() == 1 ) {
            """
            filter_tnscope_unpaired.py --vcf $vcf --out ${vcf}.vcf

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version 2>&1 | sed -e 's/Python //g')
            END_VERSIONS
            """
        }

    stub:
        """
        touch ${vcf}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process COYOTE_SEGMENTS_BED {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), val("coyote_segments"), file("*.cn-segments.raw.bed"), emit: raw_segments_bed
        path "versions.yml",                                                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        coyote_segmentator.py --vcf $vcf --sample-id ${meta.id} --segments-out ${prefix}.cn-segments.raw.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.cn-segments.raw.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process COYOTE_SEGMENTS {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(intersect_bed)

    output:
        tuple val(group), val(meta), file("*.cn-segments.panel.bed"),  emit: filtered
        tuple val(group), val(meta), file("*.cn-segments.bed"),        emit: raw
        path "versions.yml",                                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""

        """
        coyote_segmentator.py --intersect-bed $intersect_bed --sample-id ${meta.id} --panel-out ${prefix}.cn-segments.panel.bed --raw-out ${prefix}.cn-segments.bed $normal $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.cn-segments.panel.bed ${prefix}.cn-segments.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process COYOTE_SEGMENTS_JSON {
    label "process_single"
    tag "${meta.id}"

    input:
        tuple val(group), val(meta), file(bed)
    
    output:
        tuple val(group), val(meta), file("*panelmatched.json"),  emit: json_panel
        path "versions.yml",                                      emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""

        """
        cnvJSON.py --bed $bed $args --id ${meta.id} --out ${prefix}_cnvs.json --panelmatched-out ${prefix}_cnvs_panelmatched.json $normal

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def normal = meta.type.equals('normal') || meta.type.equals('N') ? "--normal" : ""
        """
        touch ${prefix}_cnvs_panelmatched.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}

process MERGE_SEGMENTS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(segments)

    output:
        tuple val(group), file("*.cn-segments.panel.merged.bed"), emit: merged

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        cat $segments > ${prefix}.cn-segments.panel.merged.bed
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}.cn-segments.panel.merged.bed
        """

}

process MERGE_JSON {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(segments)

    output:
        tuple val(group), file("*.cnvs.merged.json"), emit: merged

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        segments = segments.join(' ')
        if (meta.id.size() >= 2) {
            """
            jq $args $segments > ${group}.cnvs.merged.json
            """
        }
        else {
            """
            cat $segments > ${group}.cnvs.merged.json
            """
        }


    stub:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        segments = segments.join(' ')
        if (meta.id.size() >= 2) {
            """
            echo jq $args $segments > ${group}.cnvs.merged.json
            """
        }
        else {
            """
            echo $segments > ${group}.cnvs.merged.json
            """
        }

}

process FILTER_MANTA {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*_manta_filtered.vcf"),     emit: filtered
        tuple val(group), val(meta), file("*_manta_bnd_filtered.vcf"), emit: bnd_filtered
        path "versions.yml",                                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args   ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        filter_manta.py --vcf $vcf --filtered-out ${prefix}_manta_filtered.vcf --bnd-out ${prefix}_manta_bnd_filtered.vcf $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_manta_bnd_filtered.vcf
        touch ${prefix}_manta_filtered.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

}

process GENEFUSE_JSON_TO_VCF {
    label "process_single"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(json)

    output:
        tuple val(group), val(meta), file("*_genefuse.vcf"),    emit: genefuse_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        genefuse_json_to_vcf.py -i ${meta.id} -j $json -o ${prefix}_genefuse.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_genefuse.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process BIOMARKERS_TO_JSON {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), file(markers)

    output:
        tuple val(group), file("*.bio.json"),   emit: biomarkers_json
        path "versions.yml",                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
        tmp = (msis ?: []) + (msip ?: []) + (hrd ?: [])
        command = tmp.join(' ')

        """
        aggregate_biomarkers.py $command --out ${prefix}.bio.json --id $group

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        msis_idx = markers.findIndexOf{ it =~ 'msi_single' }
        msip_idx = markers.findIndexOf{ it =~ 'msi_paired' }
        hrd_idx = markers.findIndexOf{ it =~ 'HRD' }
        // find biomarkers //
        msis = msis_idx >= 0 ? markers[msis_idx].collect {'--msi_s ' + it} : null
        msip = msip_idx >= 0 ? markers[msip_idx].collect {'--msi_p ' + it} : null
        hrd = hrd_idx >= 0 ? markers[hrd_idx].collect {'--hrd ' + it} : null
        tmp = (msis ?: []) + (msip ?: []) + (hrd ?: [])
        command = tmp.join(' ')

        """
        echo $command
        touch ${prefix}.bio.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}

process POST_ANNOTATION_FILTERS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf)
        
    output:
        tuple val(group), val(meta), file("*.final.filtered.vcf"),    emit: filtered_vcf
        path "versions.yml",                                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix  = task.ext.prefix   ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        post_annotation_filtering.py --vcf $vcf $args --out ${prefix}.final.filtered.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
    stub:
        def prefix  = task.ext.prefix ?: "${group}"
        def args    = task.ext.args     ?: ''
        """
        touch ${prefix}.final.filtered.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
        
}
