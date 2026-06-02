# Auxiliary files

A certain number of files are required to run the pipeline, and are successively specified in the `nextflow_hopper.config`.
Only the auxiliary files needed for the profiles currently in use are listed below.

| Parameter                         | Function         | Format          | Description                                      |
|-----------------------------------|------------------|-----------------|--------------------------------------------------|
| **params/PARP_inhib/solid/GMSHem**|                  |                 |                                                  | 
| `idSnp_bed`                       |                  | `bed`           |                                                  | 
| `idSnp_bed_gz`                    |                  | `tsv`           |                                                  | 
| `idSnp_std_bed_gz`                |                  | `bed`           |                                                  | 
| **params/solid**                  |                  |                 |                                                  | 
| `mantafusions`                    | FUSIONS          | `bed`           |                                                  |
| **params/GMSHem**                 |                  |                 |                                                  | 
| `pindel_regions_bed`              | SNV-CALLING      | `bed`           |                                                  |
| **hg38 profile**                  |                  |                 |                                                  | 
| `genome_file`                     | Reference genome | `fasta`         | The main GRCh38 reference                        |
| `GENOMEDICT`                      | Reference genome | `dict`          | Dict file for the FASTA file [how was it build?] |
| `split_ref`                       | Reference genome |                 |                                                  |
| `CADD`                            | Annotation (VEP) | `tsv`           |                                                  | 
| `VEP_FASTA`                       | Annotation (VEP) | `fasta`         |                                                  | 
| `GNOMAD`                          | Annotation (VEP) | `vcf`           |                                                  | 
| `COSMIC`                          | Annotation (VEP) | `vcf`           |                                                  | 
| `gene_regions`                    | Annotation (VEP) | `bed`           |                                                  | 
| `gencode_genes`                   | Annotation (VEP) | `bed`           |                                                  | 
| `GATK_GNOMAD`                     | Annotation (VEP) |                 |                                                  | 
| `refflat`                         | Annotation (VEP) | `txt`           |                                                  |
| `blacklist`                       | Annotation (VEP) | `bed`           |                                                  | 
| `priors`                          | Annotation (VEP) |                 |                                                  | 
| `gene_gtf`                        | Annotation (VEP) | `bed`           |                                                  | 
| `verifybamidloci`                 | Annotation (VEP) | `vcf`           |                                                  | 
| **PARP_inhib**                    |                  |                 |                                                  | 
| `PON_freebayes`                   | PON              |                 |                                                  | 
| `PON_vardict`                     | PON              |                 |                                                  |
| `regions_bed`                     |                  | `bed`           |                                                  |
| `interval_list`                   |                  |                 |                                                  |
| `regions_proteincoding`           |                  | `bed`           |                                                  |
| `markgermline`                    | SNV calling      | `json`          |                                                  |
| `panel_cna`                       | CNV calling      |                 |                                                  |
| `bedgz`                           | MANTA            | `bed`           |                                                  |
| `cnvkit_reference`                | CNV calling      |                 |                                                  |
| `gatk_intervals_full`             | CNV calling      | `intervall_list`|                                                  |
| `GATK_pon`                        | CNV calling      | `hdf5`          |                                                  |
| `gatk_intervals`                  | CNV calling      | `intervall_list`|                                                  |
| **solid**                         |                  |                 |                                                  | 
| `PON_freebayes`                   | PON              |                 |                                                  |
| `PON_vardict`                     | PON              |                 |                                                  |
| `regions_bed`                     |                  | `bed`           |                                                  |
| `regions_bed_exons`               |                  | `bed`           |                                                  |
| `regions_bed_backbone`            |                  | `bed`           |                                                  |
| `regions_backbone_idsnps`         |                  | `bed`           |                                                  |
| `interval_list`                   |                  | `intervall_list`|                                                  |
| `regions_bed_qc`                  |                  | `bed`           |                                                  | 
| `interval_list_qc`                |                  | `intervall_list`|                                                  |
| `regions_proteincoding`           |                  | `bed`           |                                                  |
| `cov_probes`                      |                  | `bed`           |                                                  |
| `markgermline`                    | SNV calling      | `json`          |                                                  |
| `panel_cna`                       | CNV calling      |                 |                                                  |
| `loqusdb_vcfs`                    | CNV calling      | `vcf`           |                                                  |
| `bedgz`                           | MANTA            | `bed`           |                                                  |
| `cnvkit_reference`                | CNV calling      | `cnn`           |                                                  |
| `cnvkit_reference_exons`          | CNV calling      | `cnn`           |                                                  |
| `cnvkit_reference_backbone`       | CNV calling      | `cnn`           |                                                  |
| `gatk_intervals_full`             | CNV calling      | `intervall_list`|                                                  |
| `GATK_pon`                        | CNV calling      | `hdf5`          |                                                  |
| `gatk_intervals`                  | CNV calling      | `intervall_list`|                                                  |
| `genefuse_reference`              | FUSIONS          | `csv`           |                                                  |
| `msi_baseline`                    | BIOMARKERS       | `list`          |                                                  |
| `msi_pro_baseline`                | BIOMARKERS       | `list_baseline` |                                                  |
| **GMSHem**                        |                  |                 |                                                  | 
| `PON_freebayes`                   | PON              |                 |                                                  | 
| `PON_vardict`                     | PON              |                 |                                                  |
| `PON_tnscope`                     | PON              |                 |                                                  |
| `regions_bed`                     |                  | `bed`           |                                                  |
| `interval_list`                   |                  | `intervall_list`|                                                  |
| `regions_bed_qc`                  |                  | `bed`           |                                                  | 
| `interval_list_qc`                |                  | `intervall_list`|                                                  |
| `regions_proteincoding`           |                  | `bed`           |                                                  |
| `regions_bed_backbone`            |                  | `bed`           |                                                  |
| `regions_backbone_idsnps`         |                  | `bed`           |                                                  |
| `markgermline`                    | SNV calling      | `json`          |                                                  |
| `panel_cna`                       | CNV calling      |                 |                                                  |
| `manta`                           | MANTA            | `bed`           |                                                  |
| `cnvkit_reference`                | CNV calling      | `cnn`           |                                                  |
| `gatk_intervals_full`             | CNV calling      | `intervall_list`|                                                  |
| `GATK_pon`                        | CNV calling      | `hdf5`          |                                                  |
| `gatk_intervals`                  | CNV calling      | `intervall_list`|                                                  |