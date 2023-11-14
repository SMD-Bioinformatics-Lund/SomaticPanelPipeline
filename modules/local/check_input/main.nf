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
        check_samplesheet.py -c ${samplesheet} -o samplecheck.txt

        if [[ -e "samplecheck.txt" ]]; then
            cp ${samplesheet} "${prefix}.checked.csv"
        else
            echo "samplecheck.txt does not exist"
        fi
        """
}