# Containers 

## Open-source

Several containers are used in this pipeline and called for in the different config files found in the `configs/` folder:

- publicly available docker container:
```bash
singularity pull docker://ensemblorg/ensembl-vep:release_103
```

## Commercial

### Sentieon

A valid license for [Sentieon](https://support.sentieon.com/quick_start/index.html#document-index) is required to run the pipeline.
In the nextflow.config file, the environment variable SENTIEON_LICENSE needs to be set to the IP + port of the license server:
```bash
env 
{
  SENTIEON_LICENSE                    = '10.139.0.101:8990'
  }
```
