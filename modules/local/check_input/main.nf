process CSV_CHECK {
    label 'process_single'

    input:
        path samplesheet

    output:
        path "*.checked.csv", emit: csv

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${samplesheet.baseName}"
        """
        check_samplesheet.py -c ${samplesheet} -o samplecheck.txt --cmd-assays-json ${projectDir}/resources/cmd_assays.json

        if [[ -e "samplecheck.txt" ]]; then
            cp ${samplesheet} "${prefix}.checked.csv"
        else
            echo "samplecheck.txt does not exist"
        fi
        """
}

process GENES_ANALYZED {
	cpus 2
	tag "$group"
	time '1h'
	memory '20 GB'

	input:
		tuple val(group), val(meta)

	output:
		tuple val(group), file("${group}.genes"), emit: genes_of_interest

	script:
		def panels = (meta.diagnosis ?: '')
			.split(/\+/)
			.collect { it.trim() }
			.findAll { it }
		def panelsJson = groovy.json.JsonOutput.toJson(panels)
		"""
		jq -r --argjson panels '${panelsJson}' '
		[x
			.[]
			| select(.name as \$p | $panels | index(\$p))
			| .genes[]
		]
		| unique
		| .[]
		' ${params.all_coyote_genepanels} > ${group}.genes || touch ${group}.genes
		"""
	stub:
		"""
		touch ${group}.genes
		"""
}