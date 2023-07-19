#!/usr/bin/env nextflow

workflow CHECK_INPUT {
	take:
		csv

	main:
		//CSV_CHECK ( csv )
		csv.splitCsv( header:true, sep:',' ).set { csvmap }

		sample     = csvmap.map { create_fastq_channel(it) }

	emit:
		sample           // channel: [ val(meta), [ reads ] ]

}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
	// create meta map
	def meta = [:]
	meta.id                 = row.id
	meta.cram               = row.cram
    meta.crai               = row.crai
    meta.bai                = row.bai
	
	// add path(s) of the fastq file(s) to the meta map
	def sample = []
	sample = [ row.id, file(row.cram), file(row.crai), file(row.bai) ]

	return sample
}