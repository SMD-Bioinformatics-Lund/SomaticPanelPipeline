process GENEFUSE {
    cpus 30
    memory '70GB'
    publishDir "${params.outdir}/${params.subdir}/fusions", mode: 'copy', overwrite: true
    time '1h'
    tag "${meta.id}"
	container = "/fs1/resources/containers/genefuse-0.8.0.sif"
	errorStrategy 'ignore'

	input:
		tuple val(group), val(meta), file(r1), file(r2)

	output:
		tuple val(group), val(meta), file("${meta.id}.genefuse.json"), emit: genefuse_json
		tuple val(group), val(meta), file("${meta.id}.html"), emit: genefuse_html
	
	script:

		"""
		genefuse -t $task.cpus \\
			-r $params.genome_file \\
			-f $params.genefuse_reference \\
			-1 $r1 \\
			-2 $r2 \\
			-h ${meta.id}.html \\
			-j ${meta.id}.genefuse.json
		"""
	stub:
		"""
		touch ${meta.id}.genefuse.json ${meta.id}.html
		"""

}