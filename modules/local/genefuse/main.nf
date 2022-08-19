process GENEFUSE {
    cpus 30
    memory '70GB'
    publishDir "${params.outdir}/${params.subdir}/fusions", mode: 'copy', overwrite: true
    time '2h'
    tag "$group"
	container = "/fs1/resources/containers/genefuse-0.8.0.sif"

	input:
		tuple val(group), val(id), val(type), file(r1), file(r2)

	output:
		tuple val(group), val(id), val(type), file("${group}.genefuse.json"), emit: genefuse_json

	"""
	genefuse -t $task.cpus \\
		-r $params.genome_file \\
		-f $params.genefuse_reference \\
		-1 $r1 \\
		-2 $r2 \\
		-h ${group}.html \\
		-j ${group}.genefuse.json
	"""


}