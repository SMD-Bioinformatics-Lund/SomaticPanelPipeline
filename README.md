# README
SomaticPanelPipeline (SPP)  is a Nextflow DSL2 pipeline for somatic variant detection and annotation in targeted gene panel sequencing data from cancer patients. It is developed and maintained by the Center of Molecular Diagnostic, Lund (CMD Lund) and is designed to run on a Slurm HPC cluster using Singularity containers.
The pipeline supports multiple cancer panel assays — myeloid, lymphoid, solid tumors and PARP inhibitor response — handles both tumor-only and tumor-normal paired samples, as well as FFPE samples.
Starting from raw FASTQ files, SPP performs alignment, quality control, SNV/indel and CNV callings, fusion detection, and biomarker estimation. Results are aggregated and imported into the [Coyote](https://github.com/SMD-Bioinformatics-Lund/coyote3) variant interpretation database.

## Table of content

* [Pipeline overview](docs/pipeline_overview.png)
* [Running the pipeline](docs/running_the_pipeline.md)
* Input:
    * [containers](docs/containers.md)
* [List of used softwares](docs/list_of_used_softwares.md)
* Outputs
