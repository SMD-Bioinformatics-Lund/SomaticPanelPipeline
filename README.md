# README
SomaticPanelPipeline (SPP)  is a Nextflow DSL2 pipeline for somatic variant detection and annotation in targeted gene panel sequencing data from cancer patients. It is developed and maintained by the Center of Molecular Diagnostic, Lund (CMD Lund) and is designed to run on a Slurm HPC cluster using Singularity containers.
The pipeline supports multiple cancer panel assays — myeloid, lymphoid, solid tumors and PARP inhibitor response — handles both tumor-only and tumor-normal paired samples, as well as FFPE samples.
Starting from raw FASTQ files, SPP performs alignment, quality control, SNV/indel and CNV callings, fusion detection, and biomarker estimation. Results are aggregated and imported into the [Coyote](https://github.com/SMD-Bioinformatics-Lund/coyote3) variant interpretation database.

## Table of content

* [Pipeline overview](docs/pipeline_overview.png)
* [Running the pipeline](docs/running_the_pipeline.md)
* Input:
    * [containers](docs/containers.md)
    * [Sentieon license](docs/containers.md)
* [List of used softwares](docs/list_of_used_softwares.md)

## Pipeline overview

The SomaticPanelPipeline

A typical workflow

1-Validate input CSV and create sample metadata
2-Downsample and/or trim reads if required by the assay
3-Align reads with BWA-MEM and process UMI tags, deduplicate and apply BQSR
4-Calculate BAM QC and coverage metrics, check ID-SNPs for sample identity
5-Call SNVs and indels in parallel across bed regions using Freebayes, VarDict, TNscope and Pindel
6-Filter SNV calls against panel-of-normals (PoN), annotate with VEP and mark germline variants
7-Call copy number variants using CNVkit and GATK, detect structural variants with Manta
8-Annotate and segment CNV calls, calculate HRD and MSI biomarkers
9-Detect DNA fusions using Genefuse and Manta
10-Aggregate all results and import into the Coyote variant interpretation database