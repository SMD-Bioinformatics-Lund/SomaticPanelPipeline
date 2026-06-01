# Containers

Two containers are used in this pipeline and called for in the nextflow_hopper.config:
- base container:
```bash
SomaticPanelPipeline_2021-06-24.sif
```
This container contains all necessary tools listed [here](list_of_used_softwares.md)
- publicly available docker container:
```bash
singularity pull docker://ensemblorg/ensembl-vep:release_103
```
