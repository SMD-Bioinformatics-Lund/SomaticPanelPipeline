process COYOTE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), val(import_types), val(import_files)

    output:
        tuple val(group), file("*.coyote"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        def samples = meta instanceof List ? meta : [meta]
        def tumor_idx = samples.findIndexOf { it.type == 'tumor' || it.type == 'T' }
        tumor_idx = tumor_idx >= 0 ? tumor_idx : 0
        def normal_idx = samples.findIndexOf { it.type == 'normal' || it.type == 'N' }
        def tumor = samples[tumor_idx]
        def normal = normal_idx >= 0 ? samples[normal_idx] : null
        def process_group = samples.size() >= 2 ? "${group}p" : group
        def environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        def access_dir = [
            cnv       : 'cnv',
            transloc  : 'fusions',
            biomarkers: 'biomarkers',
            cnvprofile: 'plots',
            lowcov    : 'QC'
        ]
        def import_order = ['cnv', 'transloc', 'biomarkers', 'cnvprofile', 'lowcov']
        def seen_import_types = [] as Set
        def import_records = [import_types, import_files].transpose()
            .findAll { record -> record[0] != 'cnvprofile' || seen_import_types.add(record[0]) }
            .sort { a, b -> import_order.indexOf(a[0]) <=> import_order.indexOf(b[0]) }
        def import_args = import_records
            .findAll { record -> access_dir.containsKey(record[0]) }
            .collect { record -> "--${record[0]} /access/${params.subdir}/${access_dir[record[0]]}/${record[1]}" }

        if (tumor.purity) {
            import_args << "--purity ${tumor.purity.toFloat()}"
        }

        def command = [
            "/data/bnf/scripts/import_DSL2_to_coyote.pl",
            "--group ${params.coyote_group}",
            "--vcf /access/${params.subdir}/vcf/${vcf}",
            "--id ${process_group}",
            "--clarity-sample-id ${tumor.clarity_sample_id}",
            "--clarity_case_id ${tumor.clarity_sample_id}",
            "--clarity_control_id ${normal?.clarity_sample_id}",
            "--build 38",
            "--gens ${tumor.id}",
            "--subpanel ${tumor.diagnosis}",
            "--clarity-pool-id ${tumor.clarity_pool_id}",
            "--clarity_case_pool_id ${tumor.clarity_pool_id}",
            "--clarity_control_pool_id ${normal?.clarity_pool_id}",
            "--sample_no ${samples.size()}",
            "--case_id ${tumor.id}",
            "--control_id ${normal?.id}",
            "--profile ${environment}",
            "--assay ${params.coyote_group}",
            "--sequencing_scope panel",
            "--omics_layer DNA",
            "--sequencing_technology Illumina",
            "--pipeline ${workflow.manifest.name}",
            "--pipeline_version ${workflow.manifest.version}",
            "--case_ffpe ${tumor.ffpe ? true : false}",
            "--case_sequencing_run ${tumor.sequencing_run ?: null}",
            "--case_reads ${tumor.reads ?: null}",
            "--case_purity ${tumor.purity ? tumor.purity.toFloat() : null}",
            "--control_ffpe ${normal ? (normal.ffpe ? true : false) : null}",
            "--control_sequencing_run ${normal?.sequencing_run ?: null}",
            "--control_reads ${normal?.reads ?: null}",
            "--control_purity ${normal?.purity ? normal.purity.toFloat() : null}",
            "--paired ${samples.size() >= 2}"
        ] + import_args

        """
        cat <<-'END_COYOTE' > ${process_group}.coyote
        ${command.join(' ')}
        END_COYOTE
        """

    stub:
        def samples = meta instanceof List ? meta : [meta]
        def process_group = samples.size() >= 2 ? "${group}p" : group
        """
        touch ${process_group}.coyote
        """
}

process COYOTE_YAML {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), val(import_types), val(import_files)

    output:
        tuple val(group), file("*.coyote3.yaml"), emit: coyote_import

    when:
        task.ext.when == null || task.ext.when

    script:
        def samples = meta instanceof List ? meta : [meta]
        def tumor_idx = samples.findIndexOf { it.type == 'tumor' || it.type == 'T' }
        tumor_idx = tumor_idx >= 0 ? tumor_idx : 0
        def normal_idx = samples.findIndexOf { it.type == 'normal' || it.type == 'N' }
        def tumor = samples[tumor_idx]
        def normal = normal_idx >= 0 ? samples[normal_idx] : null
        def process_group = samples.size() >= 2 ? "${group}p" : group
        def environment = params.dev ? 'development' : params.validation ? 'validation' : params.testing ? 'testing' : 'production'
        def access_dir = [
            cnv       : 'cnv',
            transloc  : 'fusions',
            biomarkers: 'biomarkers',
            cnvprofile: 'plots',
            lowcov    : 'QC',
            cov       : 'QC'
        ]
        def import_order = ['cnv', 'transloc', 'biomarkers', 'cnvprofile', 'lowcov', 'cov']
        def yaml_lines = [
            "---",
            "subpanel: '${tumor.diagnosis}'",
            "name: '${process_group}'",
            "clarity_case_id: '${tumor.clarity_sample_id}'",
            "clarity_control_id: '${normal?.clarity_sample_id}'",
            "clarity_case_pool_id: '${tumor.clarity_pool_id}'",
            "clarity_control_pool_id: '${normal?.clarity_pool_id}'",
            "genome_build: 38",
            "vcf_files: /access/${params.subdir}/vcf/${vcf}",
            "sample_no: ${samples.size()}",
            "case_id: '${tumor.id}'",
            "control_id: '${normal?.id}'",
            "profile: '${environment}'",
            "assay: '${params.coyote_group}'",
            "sequencing_scope: 'panel'",
            "omics_layer: 'DNA'",
            "sequencing_technology: 'Illumina'",
            "pipeline: '${workflow.manifest.name}'",
            "pipeline_version: ${workflow.manifest.version}",
            "case_ffpe: ${tumor.ffpe ? true : false}",
            "case_sequencing_run: '${tumor.sequencing_run ?: null}'",
            "case_reads: ${tumor.reads ?: null}",
            "case_purity: ${tumor.purity ? tumor.purity.toFloat() : null}",
            "control_ffpe: ${normal ? (normal.ffpe ? true : false) : null}",
            "control_sequencing_run: '${normal?.sequencing_run ?: null}'",
            "control_reads: ${normal?.reads ?: null}",
            "control_purity: ${normal?.purity ? normal.purity.toFloat() : null}",
            "paired: ${samples.size() >= 2}"
        ]

        def seen_import_types = [] as Set
        def import_records = [import_types, import_files].transpose()
            .findAll { record -> record[0] != 'cnvprofile' || seen_import_types.add(record[0]) }
            .sort { a, b -> import_order.indexOf(a[0]) <=> import_order.indexOf(b[0]) }

        yaml_lines += import_records
            .findAll { record -> access_dir.containsKey(record[0]) }
            .collect { record -> "${record[0]}: /access/${params.subdir}/${access_dir[record[0]]}/${record[1]}" }

        if (tumor.purity) {
            yaml_lines << "purity: ${tumor.purity.toFloat()}"
        }

        """
cat <<'END_YAML' > ${process_group}.coyote3.yaml
${yaml_lines.join('\n')}
END_YAML
        """

    stub:
        def samples = meta instanceof List ? meta : [meta]
        def process_group = samples.size() >= 2 ? "${group}p" : group
        """
        touch ${process_group}.coyote3.yaml
        """
}
