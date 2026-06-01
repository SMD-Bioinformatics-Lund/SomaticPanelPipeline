# Getting started

Requirements for running the pipeline:
* source code
* a nextflow.config file
* Dependencies
    * [Nextflow/23.04.2](https://www.nextflow.io/)
    * [Singularity/3.8.0](https://docs.sylabs.io/guides/4.4/user-guide/#)
    * Java/13.0.2
    * perl: 5.26.2
* [auxiliary files](docs/auxiliary_files.md)

## Source code

Clone the source code from [GitHub](https://github.com/SMD-Bioinformatics-Lund/SomaticPanelPipeline#)

```bash
git clone https://github.com/SMD-Bioinformatics-Lund/SomaticPanelPipeline.git
```

## The nextflow.config file

The main config file is called `nextflow_hopper.config` and can be found in the [configs/](configs/) folder. It contains various parameters and auxiliary files required to run the pipeline. It should be copied next to `main.nf` and renamed as `nextflow.config`.