# Getting started

Requirements for running the pipeline:
* source code
* a nextflow.config file
* Dependencies
    * [Nextflow/23.04.2](https://www.nextflow.io/)
    * [Singularity/3.8.0](https://docs.sylabs.io/guides/4.4/user-guide/#)
    * Java/13.0.2
    * perl: 5.26.2 (unsure, might be in a container)
* [Auxiliary files](list currently unavailable)

## Source code

Clone the source code from [GitHub](https://github.com/SMD-Bioinformatics-Lund/SomaticPanelPipeline#)

```bash
git clone https://github.com/SMD-Bioinformatics-Lund/SomaticPanelPipeline.git
```

## The nextflow.config file

The main config file is called `nextflow_hopper.config` and can be found in the [configs/](configs/) folder. It contains various parameters and auxiliary files required to run the pipeline. It should be copied next to `main.nf` and renamed as `nextflow.config`.

## Dependencies

### Nextflow

Nextflow is a workflow orchestration framework that facilitates the creation, execution, and maintenance of reproducible computational pipelines. Its flexibility allows workflows to be deployed across different computing environments with minimal modifications.
The pipeline is implemented using the DSL2 syntax.
To run Nextflow, you will also need to have Java installed.

### Singularity

This pipeline uses Singularity containers to encapsulate software dependencies and execution environments.

### Loading dependencies on a cluster

To run the pipeline on a cluster, you will usually have to load the dependencies first as follow:
```bash
module load singularity
module load Java
module load nextflow
```

## Running the pipeline

### Running in stub-run mode

The `stub-run` (or `stub`) nextflow command allows to test the workflow logic of the pipeline, as a "dry-run", without executing the real scripts, and outputing dummy files. The stub block must be present in each nextflow process, otherwise the script section will be executed instead.

```bash
nextflow run main.nf \
        -stub-run \
        -entry SPP \
        -profile hg38,<name_of_profile_from nextflow.config> \
        --csv <path_to_input_csv_file>
```

Note that the profile `hg38` needs to always be included with the profile of analysis you want to run.

### Running the pipeline with real data
```bash
nextflow run main.nf \
        -entry SPP \
        -profile hg38,<name_of_profile_from nextflow.config> \
        --csv <path_to_input_csv_file>
```
Note that the profile hg38 needs to always be included with the profile of analysis you want to run.
Other usefull commands:
- the command `-resume`if you want to re-run the pipeline and only execute the missing/changed steps from a previous run.
- the `--dev` params allows you to redirect the results directory to results_dev