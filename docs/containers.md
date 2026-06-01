# Containers

Several containers are used in this pipeline and called for in the different config files found in the `configs/` folder:
- base container:
```bash
SomaticPanelPipeline_2021-06-24.sif
```
This container contains all necessary tools listed [here](list_of_used_softwares.md)
- publicly available docker container:
```bash
singularity pull docker://ensemblorg/ensembl-vep:release_103
```



sentieon_202112_rnaseq-expr.sif
gatk_4.1.9.0.sif
msisensor-pro-1.2.0.sif
scarHRD.sif
svdb_2.8.2.sif
cnvkit099.sif
wgs_active.sif
genefuse-0.8.0.sif
snpeff_4.3.1t2.sif
perl-json.sif
python_d4_bedtools.sif
verifybamid2_2.0.1.sif
perl-gd.sif
peddy_0.4.8.sif
python_cmdvcf0.1_pysam0.22.1.sif
seqtk_1.3.sif
fastp_0.23.sif
perl_xmljson.sif
python_cmdvcf_pysam0.22.1.sif
melt_2.2.2.sif
entieon_20250520.sif => DNAscope
multiqc_1_14_0.img