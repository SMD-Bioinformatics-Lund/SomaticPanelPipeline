process CSV_CHECK {

    input:
        path samplesheet

    output:
        path '*.csv'       , emit: cool


    script: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
        """
        echo bla > bla.csv
        """

    stub: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
        """
        echo bla > bla.csv
        """
}