process CSV_CHECK {

    input:
    path samplesheet

    output:
    path "${samplesheet.baseName}.checked.csv", emit: csv

    script:
    """
    check_samplesheet.py -c ${samplesheet} -o samplecheck.txt

    if [[ -e "samplecheck.txt" ]]; then
        cp ${samplesheet} "${samplesheet.baseName}.checked.csv"
    else
        echo "samplecheck.txt does not exist"
    fi
    """

    stub:
     """
		echo ${samplesheet} 
        touch ${samplesheet}.checked.csv
     """
}