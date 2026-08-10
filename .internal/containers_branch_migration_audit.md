# Containers Branch Migration Audit

Date: 2026-08-09

Scope: static comparison of the current working tree on branch `containers` against `main` for `/data/bnf/dev/ram/Pipelines/DSL2/SomaticPanelPipeline`.

I did not run Nextflow, execute containers, copy files to `/fs1`, or modify pipeline code. This report is based on file diffs, process signatures, config selectors, includes, output declarations, and publish rules.

## Executive Summary

The migration is substantial and is not yet safe to call result-equivalent to `main` from static review alone.

The broad module split is mostly wired into workflows and configs: the new CNVkit and GATK subworkflow config files are included by `configs/nextflow.hopper.config`, many removed monolithic process names have matching split process selectors, and several tool-specific modules have been introduced.

However, exact output equivalence is not guaranteed. I found multiple high-risk areas where output channel shapes, published files, or process semantics changed. Some changes are expected from the refactor, but several need deliberate confirmation or cleanup before claiming the same results.

## Branch State

Current branch: `containers`

The working tree is dirty relative to `main`. The audit compared `main` to the current working tree, not only `HEAD`.

Scale of changes relative to `main`:

- 98 files changed
- About 7,098 insertions and 2,461 deletions
- Large changes in `bin/`, `configs/modules/`, `modules/local/`, `subworkflows/local/`, and `workflows/`

## Config Migration

### Looks Migrated

The new split config files are included from `configs/nextflow.hopper.config`:

- `configs/modules/cnvkit_cnv_calling.config`
- `configs/modules/gatk_cnv_calling.config`
- Existing module configs remain included, including `cnv_calling.config`, `snv_calling.config`, `snv_annotate.config`, `qc.config`, `idsnp.config`, and others.

The following split selector namespaces are present and broadly match the new process tree:

- `.*CNV_CALLING:CNVKIT:*`
- `.*CNV_CALLING:GATK:*`
- `.*SNV_CALLING:*`
- `.*SNV_ANNOTATE:*`
- `.*QC:*`
- `.*ID_SNP:*`

### Config Risks

1. `configs/modules/cnv_calling.config` now only covers Manta and merge/filter steps around Manta. CNVkit and GATK settings moved out correctly, but this means future edits must happen in three files:
   - `cnv_calling.config`
   - `cnvkit_cnv_calling.config`
   - `gatk_cnv_calling.config`

2. Some process containers changed from old broad containers to small tool containers. That is expected, but exact result equivalence depends on tool version parity and Python script parity. Static review cannot prove identical outputs.

3. Some config blocks publish files that were not published before, or publish through a different process than before. See "Published Output Capture" below.

## Workflow And Subworkflow Migration

### Looks Migrated

The following major splits are wired into workflows/subworkflows:

- `CNV_CALLING` now delegates to nested `CNVKIT` and `GATK` subworkflows.
- `SNV_CALLING` now uses tool-specific modules for `vcflib`, `bcftools`, `vcftools`, `tabix`, `bedtools`, and caller-specific filters.
- `BAM_QC` now uses separate `sambamba`, `bcftools`, `bedtools`, `idSnp`, and `qc` processes.
- `SNV_ANNOTATE` now uses separate `vep`, `vcfanno`, `bedtools`, `tabix`, and filter modules.

### High-Risk Workflow Changes

1. `SNV_CALLING` channel shape changed widely.

   Examples:

   - `FREEBAYES` old output:
     `val("freebayes"), val(group), file("freebayes_${bed}.vcf")`

   - `FREEBAYES` current output:
     `val(group), val(meta), val("freebayes"), file("freebayes_${bed}.vcf.raw")`

   That is a good direction for metadata preservation, but it means all downstream `groupTuple`, `join`, and file-name logic must be checked carefully. Static review shows the downstream code was updated, but exact equivalence is not proven.

2. `SNV_CALLING` final aggregation order has changed.

   The current code sorts VCF parts by file name before concatenation. That is likely needed for determinism, but if `main` relied on original channel emission order, final VCF ordering may differ.

3. `SNV_ANNOTATE` CNV backbone filter currently passes `params.regions_backbone_idsnps` directly into a process input declared as `path(intersect)`.

   The empty-value guard is present for `gmshem` and `solid`, which avoids the earlier empty `file()` error. But because the process input is `path(intersect)`, this should be verified for staging behavior. The safer pattern is to only build a file/path channel inside the profile condition.

4. `BAM_QC` has been split and now routes low coverage through:

   `SAMBAMBA_PANEL_DEPTH -> PANEL_DEPTH -> LOWCOV_INTERSECT -> LOWCOV`

   This matches the tool split direction, but the LOWCOV input changed from BAM/BAI in `main` to an intersected BED in the current branch. Output equivalence depends on `panel_depth.py`, `bedtools intersect`, and `overlapping_genes.py` producing the same `.lowcov.bed` as the old combined implementation.

5. `ID_SNP` handling changed from one `ALLELE_CALL` process to:

   `BCFTOOLS_MPILEUP_CALL -> BCFTOOLS_ANNOTATE_IDSNP -> BCFTOOLS_QUERY_IDSNP -> GENOTYPE2JSON -> SNP_CHECK -> PAIRGEN_CDM`

   This is structurally cleaner, but `SNP_CHECK` output now carries `meta`, and the JSON naming depends on `idsnp_controller.py`. This area already produced runtime failures earlier and should be treated as high risk until tested.

## Process Signature Changes

These process signatures changed in ways that affect downstream channels:

- `modules/local/sentieon/main.nf`
  - `BWA_UMI` no longer emits `bam_umi_markdup`; split processes now provide sort/locus/markdup.
  - `MARKDUP` now requires score and score index inputs.
  - `BQSR_UMI` output emit changed from `bam_varcall` to `bam_bqsr`.
  - `TNSCOPE` now emits caller label `tnscope`.

- `modules/local/cnvkit/main.nf`
  - `CNVKIT_CALL` no longer emits VCF or nexus/BAF-logR outputs; these moved to `CNVKIT_EXPORT_VCF` and `CNVKIT_EXPORT_NEXUS_OGT`.
  - `CNVKIT_GENS` now emits uncompressed `.bed`; compression/indexing moved to tabix processes.

- `modules/local/GATK/main.nf`
  - Old `GATKCOV_COUNT`, `GATKCOV_CALL`, and `FILTER_MERGE_GATK` were split into many GATK command processes.
  - `GATK2VCF` input order changed from `[group, seg, meta]` to `[group, meta, seg]`.
  - `POSTPROCESS` now consumes extracted shard paths instead of tar inputs.

- `modules/local/qc/main.nf`
  - `CONTAMINATION` no longer emits `.contaminationpy`; that moved to `CONTAMINATION_CDM`.
  - `CONTAMINATION` now emits `.value` and keeps `meta`.
  - `LOWCOV` input changed from BAM/BAI to BED.
  - `LOWCOV_D4` now emits `meta` in output tuples.

- `modules/local/idSnp/main.nf`
  - `SNP_CHECK` and `PAIRGEN_CDM` now carry `meta` in outputs.

These changes can be valid, but they are not result-neutral unless every downstream consumer was updated and old published files are still produced with the same names.

## Published Output Capture

### Likely Preserved

The main user-facing output folders are still represented:

- `${outdir}/${subdir}/vcf`
- `${outdir}/${subdir}/svvcf`
- `${outdir}/${subdir}/svvcf/merged`
- `${outdir}/${subdir}/cnvkit`
- `${outdir}/${subdir}/cnvkit/segments`
- `${outdir}/${subdir}/gatkcov`
- `${outdir}/${subdir}/plots`
- `${outdir}/${subdir}/QC`
- `${outdir}/${subdir}/QC/contamination`
- `${outdir}/${subdir}/cov`
- `${outdir}/cron/gens`
- `${crondir}/qc`
- `${crondir}/idsnp`
- `${crondir}/contamination`

### High-Risk Publish Differences

1. CNVkit `.cnvkit` files are now explicitly published.

   `configs/modules/cnvkit_cnv_calling.config` publishes `*.cnvkit` from `CNVKIT_EXPORT_NEXUS_OGT`.

   In `main`, `CNVKIT_CALL` published `*call*.cns` and `*.vcf`, but not obviously `*.cnvkit`. This may be an added published output, not a preserved one.

2. CNVkit GENS publishing moved.

   In `main`, `CNVKIT_GENS` published `*.bed.gz` directly to `${outdir}/${subdir}/gens`.

   Current branch:
   - `CNVKIT_GENS` creates `.bed`
   - tabix/bgzip processes compress/index
   - `MERGE_GENS` publishes `*.bed.gz*`

   This can work because inputs are staged into `MERGE_GENS`, but it means the user-facing `.bed.gz*` capture depends on `MERGE_GENS` running and re-emitting staged inputs. That is a behavioral change in where publishing happens.

3. SNV final normalized VCF publishing moved from old `CONCATENATE_VCFS`/`AGGREGATE_VCFS` style to `BCFTOOLS_NORM_EXACT` and `VCF_AGGREGATE_SORT`.

   The final `vcf` folder still has publish rules, but exact filenames and intermediate capture should be checked. Current `BCFTOOLS_NORM` naming prepends caller names and uses `_norm.vcf.gz`; the second norm also uses `_norm.vcf.gz`. This may not match old final file names unless config prefixes are exactly compensating.

4. `CONTAMINATION` and `CONTAMINATION_CDM` publishing is split.

   Main published PNG/TXT and `.contaminationpy` from `CONTAMINATION`.

   Current branch publishes PNG/TXT from `CONTAMINATION` and `.contaminationpy` from `CONTAMINATION_CDM`. That is structurally acceptable, but only if both always run for all cases where the old process created all files.

5. Some intermediate process outputs are now not published where the old monolithic process published all produced files.

   This is probably desired for smaller modules, but it means "published files exactly same as before" must be checked file-by-file. Static review cannot confirm this across all profiles.

## Exact Result Equivalence

I cannot say the current branch will produce exactly the same results as `main`.

Reasons:

1. Many Perl or shell paths were replaced by Python scripts. Even if intended to mimic behavior, exact equivalence requires fixtures or golden-output comparison.

2. Some tool versions changed with smaller containers:
   - FreeBayes/vcflib/bcftools are now explicit newer tool containers in some places.
   - CNVkit container changed from `cnvkit099.sif` to `cnvkit-0.9.9.sif`.
   - Python utility containers replaced broad historical pipeline containers.

3. File sorting and channel grouping changed. This may alter record ordering in concatenated/aggregated VCFs.

4. Several output channel shapes changed by adding `meta`, `vc`, or `part`. That is needed for the refactor, but downstream joins can silently change multiplicity if keys are not exact.

5. Multiple previously monolithic processes are now split. Any missing command, missing option, or altered filename between split stages can change results.

## Specific Findings To Fix Or Verify

### Must Verify Before Claiming Equivalence

1. `SNV_CALLING` final VCF names and contents.

   Check that these old outputs still exist with the same names:

   - caller-concatenated VCFs
   - normalized caller VCFs
   - aggregate VCF
   - sorted aggregate VCF

2. `CNVKIT` published output set.

   Confirm whether `*.cnvkit` should be published. If not published in `main`, remove or restrict that publish rule.

3. `CNVKIT_GENS` `.bed.gz*` publishing.

   Current publishing through `MERGE_GENS` is indirect. Confirm all old `.baf.bed.gz`, `.cov.bed.gz`, `.tbi`, `.gens`, and `.gens_v4_somatic` files are captured with the same names.

4. `BAM_QC` lowcov output.

   Confirm new `sambamba depth -> panel_depth.py -> bedtools intersect -> LOWCOV` produces the same `.lowcov.bed` as the old LOWCOV process.

5. `ID_SNP` outputs.

   Confirm `*.genotypes.json`, `*.T.idsnp`, `*.N.idsnp`, and `*.pairgen` are created with the same names and content as before.

6. `CONTAMINATION` outputs.

   Confirm Python replacement creates the same `.txt`, `.png`, `.value`, and `.contaminationpy` behavior as the old Perl process for tumor-only and tumor-normal inputs.

7. `SNV_ANNOTATE` CNV backbone input staging.

   `CNV_BACKBONE_FILTER` currently receives `params.regions_backbone_idsnps` directly while the module expects `path(intersect)`. Verify this stages the bed file correctly. If not, build a path channel only inside the `gmshem/solid` branch.

### Lower-Risk But Worth Cleaning

1. Remove commented obsolete includes and old aggregation comments in `SNV_CALLING` once the new path is accepted.

2. Normalize process naming and config selector comments where old process names remain in comments.

3. Review whether all new Python scripts have shebangs if they are called directly as executables.

4. Make sure all optional outputs are intentional. Optional outputs can hide missing expected files if used too broadly.

## Deep Static Audit Addendum

Date: 2026-08-10

This addendum narrows the audit to command preservation, published-file preservation, output-name collisions, and process-by-process risks. It is still static review only: I did not run Nextflow, containers, or test data.

## Published Output Differences Against `main`

### Newly Published Or More Broadly Published

These are files that are now explicitly captured by a publish rule but were not clearly published by the old monolithic process/config path in `main`.

| Area | Current process/config | Current published files | Old baseline | Assessment |
|---|---|---|---|---|
| CNVkit nexus export | `CNV_CALLING:CNVKIT:CNVKIT_EXPORT_NEXUS_OGT` in `configs/modules/cnvkit_cnv_calling.config` | `*.cnvkit` to `${outdir}/${subdir}/cnvkit` | Old `CNVKIT_CALL` publish patterns were `*call*.cns` and `*.vcf`; `*.cnvkit` was not obviously published | Newly published unless old output declaration/publish caught it indirectly. If commercial output contract should be unchanged, confirm whether this file belongs in published output. |
| GATK collected counts | `CNV_CALLING:GATK:GATK_COLLECT_READ_COUNTS` | `*.hdf5` to `${outdir}/${subdir}/gatkcov` | Old `GATKCOV_COUNT` output/publish focused on standardized and denoised copy-ratio TSVs | Likely newly published intermediate. If only final user-facing files should be published, remove this publish rule or narrow it. |
| GATK denoised plot | `CNV_CALLING:GATK:GATK_PLOT_DENOISED_COPY_RATIOS` | `*.denoised.png` to `${outdir}/${subdir}/plots` | Old `GATKCOV_COUNT` command created the plot, but it was not a declared output in the old process | Newly published useful QC plot, but not baseline-equivalent. |
| GATK model-segment intermediates | `CNV_CALLING:GATK:GATK_MODEL_SEGMENTS` | `*.cr.seg`, `*.hets.tsv`, `*.modelFinal.seg` to `${outdir}/${subdir}/gatkcov` via `*.{seg,tsv}` | Old `GATKCOV_CALL` published `*.seg`, so `.cr.seg` and `.modelFinal.seg` may have been captured if declared; old process output declaration should be checked. `*.hets.tsv` is very likely newly published | Publish rule is too broad for strict parity. Split it into old-public outputs and internal-only outputs. |
| CNVkit GENS indexes | `CNV_CALLING:CNVKIT:MERGE_GENS` | `*.bed.gz*` to `${outdir}/${subdir}/gens` | Old `CNVKIT_GENS` published `*.bed.gz`; old `MERGE_GENS` published `*.bed.gz*` | If staged per-part `.bed.gz.tbi` files are emitted by `MERGE_GENS`, indexes may now be published more broadly than before. |
| Split VCF utility outputs | `SNV_CALLING:BCFTOOLS_NORM_EXACT`, `VCF_AGGREGATE_SORT`, `FILTER_PINDEL` | Final normalized and filtered VCFs to `${outdir}/${subdir}/vcf` | Old publish was from `CONCATENATE_VCFS`, `AGGREGATE_VCFS`, and `PINDEL_CALL` | Final files are probably preserved, but intermediate published set and exact names changed. See below. |

### No Longer Published Or At Risk

These were published by old `main` process selectors and are now either unpublished, renamed, or only indirectly published by a later split process.

| Area | Old process/config | Old published files | Current state | Assessment |
|---|---|---|---|---|
| SNV caller concatenation | `SNV_CALLING:CONCATENATE_VCFS` | All declared non-`versions.yml` outputs to `${outdir}/${subdir}/vcf` | Current selector is commented out; split path uses `VCF_CONCAT`, `VCF_CONCATENATED_SORT`, `TABIX_BGZIPTABIX`, `BCFTOOLS_NORM`, `BCFTOOLS_NORM_EXACT` | Intermediate concat/sort files from the old monolithic process are no longer published unless recreated by later steps. Confirm which final file names are contractually required. |
| SNV aggregate | `SNV_CALLING:AGGREGATE_VCFS` | Aggregate VCF outputs to `${outdir}/${subdir}/vcf` | Current `BT_AGG` is internal; `VCF_AGGREGATE_SORT` publishes `${group}.agg...` outputs | Likely preserved final sorted aggregate, but unsorted aggregate may no longer be published. |
| Raw Pindel files | `SNV_CALLING:PINDEL_CALL` | All declared non-version outputs to `${outdir}/${subdir}/vcf` | Current `PINDEL_CALL` and `PINDEL_TO_VCF` have no publish rule; only `FILTER_PINDEL` publishes `*_pindel.vcf` | Raw Pindel intermediate outputs are no longer published. This is probably desirable, but not baseline-equivalent. |
| CNVkit per-part GENS output | `CNV_CALLING:CNVKIT_GENS` | `*.bed.gz` directly to `${outdir}/${subdir}/gens` | Current `CNVKIT_GENS` emits uncompressed `.bed`; bgzip/tabix split; `MERGE_GENS` publishes `*.bed.gz*` | Final compressed files can be preserved, but publishing now depends on `MERGE_GENS` receiving/staging those files. If `MERGE_GENS` is skipped, the old GENS files are not published. |
| Contamination CDM file | `QC:CONTAMINATION` | `*.contaminationpy` to `${crondir}/contamination` | Current `QC:CONTAMINATION` publishes only `*.png` and `*.txt`; `QC:CONTAMINATION_CDM` publishes `*.contaminationpy` | Final file can be preserved if `CONTAMINATION_CDM` always runs after `CONTAMINATION`. This is a split-process dependency now. |
| ID-SNP final VCF in QC namespace | `QC:ALLELE_CALL` | In old QC config only `*.genotypes.json` was published; old ID_SNP config published `.final.vcf` under `ID_SNP:ALLELE_CALL` | Current `ID_SNP:BCFTOOLS_ANNOTATE_IDSNP` publishes `.final.vcf`, but `QC:BCFTOOLS_ANNOTATE_IDSNP` does not | If `BAM_QC` is the active namespace, `.final.vcf` is not published from QC unless the QC config is updated. If only `ID_SNP` workflow needs it, current idsnp config is OK. |
| CNVkit `.vcf` from call process | `CNV_CALLING:CNVKIT_CALL` | `*.vcf` to `${outdir}/${subdir}/svvcf` | Current moved to `CNVKIT_EXPORT_VCF` and `CNVKIT_EXPORT_VCF_TC` | Good split if names match old. Must verify TC and non-TC output names do not diverge. |

## Output Name Collision Risks

| Risk | Current code | Why it can collide | Correction |
|---|---|---|---|
| CNVkit TC/non-TC call names | `modules/local/cnvkit/main.nf:CNVKIT_CALL` writes `${prefix}.${part}.call.cns` unless `meta.purity && tc == "true"` | If `CNVKIT_CALL_TC` runs with `tc == "true"` but no usable `meta.purity`, it produces the same `${prefix}.${part}.call.cns` as the non-TC call. Both publish to `cnvkit/segments` with `*call*.cns` | Only run TC branch when purity is present, or force TC output to use a distinct suffix from config, e.g. `.call.tc.cns` even when purity is absent. |
| Generic BEDTOOLS output | The same `BEDTOOLS_INTERSECT` process is included under several workflow names, for example `SNV_CALLING:BEDTOOLS_INTERSECT`, `SNV_ANNOTATE:CNV_BACKBONE_FILTER`, `SNV_ANNOTATE:FILTER_FOR_CNV`, `CNV_ANNOTATE:GENE_INTERSECT`, and `CNV_ANNOTATE:COYOTE_SEGMENTS_INTERSECT`. Each configured workflow name sets `task.ext.prefix` and `task.ext.suffix`. | This is not an active collision in the checked configs, because the current workflow names set distinct prefixes/suffixes. It becomes a risk only if a future include of the same bedtools process is added without a unique `ext.prefix`/`ext.suffix`, because the module default is `${meta.id}.bed`. | Keep the current explicit `ext.prefix`/`ext.suffix` selectors for every workflow name that uses `BEDTOOLS_INTERSECT`. When adding a new `include { BEDTOOLS_INTERSECT as ... }`, add a matching `withName` block for that workflow process name. |
| `BEDTOOLS_INTERSECT` dual emits | Same process emits `intersected` and `vcf_intersected` with the same physical file but different tuple shapes. | This is not a filename collision. It is a channel-contract risk: one emit carries `[group, meta, file]`, the other carries `[group, meta, vc, file]`. A downstream consumer can accidentally use the wrong emit and lose or add the `vc` caller label. | Keep both emits only while both tuple shapes are genuinely needed. Document the expected emit in each subworkflow, or standardize the tuple shape after the current refactor settles. |
| `MARK_GERMLINES` paired prefix | `mark_germlines.py ... --out ${prefix}p.agg.pon.vep.markgerm.vcf` for paired samples | If `group` or `prefix` already carries a paired suffix, output can become `...pp...`. It also differs from tumor-only `${prefix}.agg...` | Move paired suffix decision into config or create a helper channel field for case id. Avoid hidden string mutation inside the module. |
| Manta filter names | `FILTER_MANTA` calls `filter_manta.py --vcf $vcf $args`; `FILTER_MANTA_TUMOR` and `FILTER_MANTA_NORMAL` now pass both `--id ${meta.id}` and `--out-prefix ${meta.id}`. | Fixed. Without `--out-prefix`, the script could still infer the same name, but explicit output prefix is clearer and prevents future drift. | No remaining filename collision expected unless duplicate sample IDs are sent to the same publish folder. |
| `TABIX_BGZIPTABIX` generic default | The generic tabix module writes `${prefix}.${vcf.getExtension()}.gz`, with `prefix = task.ext.prefix ?: "${vcf.baseName}"`. Current uses include `SNV_CALLING:TABIX_BGZIPTABIX`, `SNV_ANNOTATE:FILTER_FOR_CNV_TABIX`, and CNVkit BAF/COV compression. | This is not an active collision when each input file basename is unique. It is a naming risk if a future workflow passes two staged files with the same basename, or passes an already compressed file and unintentionally creates `.gz.gz`. | Use the input basename default only where the incoming file names are already the desired output names. Add `ext.prefix` in the workflow-specific config selector when the compressed output name is part of the public contract. Stub naming has been aligned with runtime naming. |
| `BCFTOOLS_NORM` naming | Runtime name uses `vc ? "${vc}_${prefix}" : "${prefix}"`. | Fixed. The stub previously used `${prefix}_${vc}`, which did not match runtime output names. | No production collision identified; this was a stub/dry-run mismatch. |
| `COYOTE_SEGMENTS_JSON` output pattern | Output catches `*panelmatched.json`; script writes names from `cnvJSON.py` using `--id`, not `ext.prefix` | Config `ext.prefix` is ignored, so changing prefix in config will not change output names | Add explicit output args to `cnvJSON.py` or pass/use `--out-prefix`. |

## Process-By-Process Command Preservation Notes

### SNV Calling

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `SNV_CALLING:FREEBAYES` | Runs FreeBayes only | Old `FREEBAYES` ran FreeBayes plus post filters and AD removal | Internal only | Good split by tool. Final equivalence depends on the follow-up `vcffilter`, `vcfglxgt`, Python filter, and `bcftools annotate` sequence. |
| `SNV_CALLING:FREEBAYES_FILTER_LOWCOV` | `vcffilter` LowCov | Part of old `FREEBAYES` | Internal only | Uses old `ext.args2` thresholds. Output name `*.lowcov.vcf` is staged in `out/`, avoiding input capture. |
| `SNV_CALLING:FREEBAYES_FILTER_LOWFRQ` | `vcffilter` LowFrq | Part of old `FREEBAYES` | Internal only | Uses old `ext.args3`. |
| `SNV_CALLING:FREEBAYES_VCFGLXGT` | `vcfglxgt` | Part of old `FREEBAYES` | Internal only | Generic vcflib process; good tool placement. |
| `SNV_CALLING:FILTER_FREEBAYES` | Python paired/unpaired FreeBayes filter | Old Perl FreeBayes filters | Internal only | Logic must match Perl. See Perl parity report. |
| `SNV_CALLING:FREEBAYES_REMOVE_AD` | `bcftools annotate -x FORMAT/AD` | Part of old `FREEBAYES` | Internal only | Correct tool placement. Name in config still says FreeBayes but process should be a generic `BCFTOOLS:ANNOTATE` alias if not already. |
| `SNV_CALLING:VARDICT` | VarDict caller | Same caller | Internal only | Filter moved to Python `FILTER_VARDICT`. |
| `SNV_CALLING:FILTER_VARDICT` | Python paired/unpaired VarDict filter | Old Perl VarDict filters | Internal only | Output naming derives from raw VCF. |
| `SNV_CALLING:TNSCOPE` | Sentieon TNscope caller | Same caller | Internal only | Still contains tumor/normal indexing inside process. Functionally OK, but not yet structurally simplified. |
| `SNV_CALLING:FILTER_TNSCOPE` | Python paired/unpaired TNscope filter | Old Perl TNscope filters | Internal only | Output pattern is broad `*.vcf`; because task workdir stages input VCF, this can capture input as well if not isolated. Needs explicit output naming or `out/` directory. |
| `SNV_CALLING:BEDTOOLS_INTERSECT` | Restricts caller VCFs to regions | Same bedtools command | Published to `vcf` | Config sets `${meta.id}_${vc}_intersected.vcf`; collision risk low if `vc` is always unique. |
| `SNV_CALLING:SVDB_MERGE_SINGLES` | Merges MELT singles | Same SVDB operation | Internal only | Only relevant when normal sample exists. |
| `SNV_CALLING:VCF_CONCAT` | Concatenates parts per caller | Part of old `CONCATENATE_VCFS` | Internal only | Old monolithic publish is removed. |
| `SNV_CALLING:VCF_CONCATENATED_SORT` | Sorts concatenated VCF | Part of old `CONCATENATE_VCFS` | Internal only | Sorting may alter record order compared with channel-emission order. |
| `SNV_CALLING:TABIX_BGZIPTABIX` | Compress/index | Part of old concat/norm path | Internal only | Generic module, but output prefix should be controlled explicitly where names matter. |
| `SNV_CALLING:BCFTOOLS_NORM` | First normalization | Part of old concat path | Internal only | Uses `${group}.first`; final name differs from old intermediates. |
| `SNV_CALLING:BCFTOOLS_NORM_EXACT` | Duplicate exact removal and index | Part of old concat path | Published to `vcf` | This appears to be the current final caller-level published VCF. Confirm filenames match old contract. |
| `SNV_CALLING:FILTER_PINDEL` | Python Pindel filter | Old Pindel filter in `PINDEL_CALL` | Published to `vcf` | Raw Pindel outputs are no longer published. |
| `SNV_CALLING:BT_AGG` | Python aggregate VCF | Old `AGGREGATE_VCFS` Perl aggregate | Internal only | Final sorted aggregate is published by `VCF_AGGREGATE_SORT`. |
| `SNV_CALLING:VCF_AGGREGATE_SORT` | Sort aggregate VCF | Part of old aggregate flow | Published to `vcf` | Final output likely preserved, but unsorted aggregate may not be. |

### SNV Annotation

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `SNV_ANNOTATE:PON_FILTER` | Python PON filter | Old Perl `filter_with_pon.pl` | Published to `vcf` | Must compare filter header/order and caller-specific INFO behavior. |
| `SNV_ANNOTATE:FFPE_PON_FILTER` | Python FFPE PON filter | Old Perl `filter_with_ffpe_pon.pl` | Published to `vcf` | Solid-only. |
| `SNV_ANNOTATE:ANNOTATE_VEP` | VEP annotation | Same VEP command | Published to `vcf` | Command args appear preserved. |
| `SNV_ANNOTATE:FIX_VEP_GNOMAD` | Python VEP gnomAD normalization | Old Perl helper likely inline or separate | Internal only | Good to keep internal if old final only consumed fixed VEP. |
| `SNV_ANNOTATE:MARK_GERMLINES` | Python germline marking | Old Perl `mark_germlines.pl` | Published to `vcf` | Paired filename suffix `p` should be verified. |
| `SNV_ANNOTATE:GERMLINE_FOR_CNVKIT` | Python germline VCF for CNVkit | Old Perl `germline_for_cnvkit.pl` | Internal only | Feeds CNVkit nexus/gens. |
| `SNV_ANNOTATE:FILTER_FOR_CNV` | Bedtools region filter on germline VCF | Old bedtools in filters module | Published to `vcf` | Prefix is hard-coded to `${group}_vardict.germlines`; this is caller-specific despite being a generic bedtools process. |
| `SNV_ANNOTATE:FILTER_FOR_CNV_TABIX` | Compress/index germline VCF | Previously inside broader command flow | Internal only | If old `.tbi` was published, current is not. |
| `SNV_ANNOTATE:CNV_BACKBONE_FILTER` | Bedtools `-v` backbone removal | Old filter module bedtools | Published to `vcf` | Passing raw params string to a `path` input remains a staging risk. |
| `SNV_ANNOTATE:POST_ANNOTATION_FILTERS` | Python/pysam final filtering | Existing Python | Published to `vcf` | Now uses `--out`, which is the right direction. |

### CNVkit CNV Calling

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `CNV_CALLING:CNVKIT:CNVKIT_FULL/BACKBONE/EXONS` | CNVkit batch by region part | Old `CNVKIT_BATCH`, `CNVKIT_BACKBONE`, `CNVKIT_EXONS` | Published to `cnvkit` | Names should preserve `${group}.${meta.id}.${part}.{cnr,cns}`. |
| `CNV_CALLING:CNVKIT:CNVKIT_PLOT` | CNVkit scatter plot | Old `CNVKIT_PLOT` | Published to `plots` | Uses germline VCF. Earlier runtime showed VCF staging/naming issue; input channel should be audited. |
| `CNV_CALLING:CNVKIT:CNVKIT_CALL` | `cnvkit.py call` | First half of old `CNVKIT_CALL` | Publishes `*call*.cns` | VCF export moved out. TC/non-TC collision risk exists. |
| `CNV_CALLING:CNVKIT:CNVKIT_EXPORT_VCF` | `cnvkit.py export vcf` | Second half of old `CNVKIT_CALL` | Publishes `*.vcf` to `svvcf` | Correct split if output names match old. |
| `CNV_CALLING:CNVKIT:CNVKIT_EXPORT_NEXUS_OGT` | `cnvkit.py export nexus-ogt` | Old command may have existed but was not clearly published | Publishes `*.cnvkit` | Newly published unless intentionally added. |
| `CNV_CALLING:CNVKIT:CNVKIT_GENS` | Python conversion from CNVkit CNR and VCF to BAF/COV BED | Old Perl `generate_gens_data_from_cnvkit.pl` plus compression | Internal only | Compression moved to tabix processes. |
| `CNV_CALLING:CNVKIT:CNVKIT_BAF_BGZIPTABIX/CNVKIT_COV_BGZIPTABIX` | Compress/index BAF/COV BED | Previously done by Perl/shell in GENS flow | Internal only | Good tool split. |
| `CNV_CALLING:CNVKIT:MERGE_GENS` | Shell helper creates merged GENS files and cron commands | Old `MERGE_GENS` | Publishes `*.bed.gz*`, `*.gens`, `*.gens_v4_somatic` | Container is set to bedtools but script also uses bgzip/tabix according to versions. Make sure container has all tools. |

### GATK CNV Calling

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `GATK:GATKCOV_BAF` | `CollectAllelicCounts` | Same | Internal only | No publish in old/current. |
| `GATK:GATK_COLLECT_READ_COUNTS` | `CollectReadCounts` | First command in old `GATKCOV_COUNT` | Publishes `*.hdf5` | Likely newly published intermediate. |
| `GATK:GATK_DENOISE_READ_COUNTS` | `DenoiseReadCounts` | Second command in old `GATKCOV_COUNT` | Publishes `*.standardizedCR.tsv`, `*.denoisedCR.tsv` | Preserves old published files. |
| `GATK:GATK_PLOT_DENOISED_COPY_RATIOS` | `PlotDenoisedCopyRatios` | Third command in old `GATKCOV_COUNT` | Publishes `*.denoised.png` | Newly published if old process did not declare it. |
| `GATK:GATK_MODEL_SEGMENTS` | `ModelSegments` | First command in old `GATKCOV_CALL` | Publishes `*.cr.seg`, `*.hets.tsv`, `*.modelFinal.seg` | Publish too broad for parity. |
| `GATK:GATK_CALL_COPY_RATIO_SEGMENTS` | `CallCopyRatioSegments` | Second command in old `GATKCOV_CALL` | Publishes `*.called.seg` | Preserved. |
| `GATK:GATK_PLOT_MODELED_SEGMENTS` | `PlotModeledSegments` | Third command in old `GATKCOV_CALL` | Publishes `*.modeled.png` | Preserved. |
| `GATK:GATK2VCF` | Python GATK seg to VCF | Old `gatk_to_vcf.py`/Perl flow | Publishes to `svvcf` | Preserve if output names match. |
| `GATK:MERGE_GATK_TUMOR` | Python tumor VCF merge | Old Perl `mergeGATK_tumor.pl` | Publishes to `svvcf` | Earlier direct execution failed due script shebang/container invocation in copied test tree; current repo script has shebang. |
| `GATK:GATK_COUNT_GERMLINE` | Germline `CollectReadCounts` | Same | Internal only | Prefix controlled by sample id. |
| `GATK:GATK_CALL_PLOIDY` | `DetermineGermlineContigPloidy` + tar | Same old process command | Internal only | Still contains `tar`, acceptable because it packages tool output. |
| `GATK:GATK_CALL_GERMLINE_CNV` | `GermlineCNVCaller` + tar | Same old process command | Publishes `*.vcf` according to config but process output is `*.tar` | This publish rule probably does nothing for this process. Final VCFs are emitted by `POSTPROCESS`. |
| `GATK:GATK_EXTRACT_GERMLINE_CNV_SHARDS` | Extract tar shards | New helper | Internal only | Split helper; not user-facing. |
| `GATK:POSTPROCESS` | `PostprocessGermlineCNVCalls` | Same old process command | Publishes three `*_gatk_*.vcf.gz` files | Preserved. |
| `GATK:FILTER_GATK` | Python GATK filter | Old Perl `filter_gatk.pl` | Publishes to `svvcf` | Preserved if Python parity holds. |
| `GATK:MERGE_GATK` | Python merge normal GATK VCF | Old Perl `mergeGATK.pl` | Publishes to `svvcf` | Preserved if Python parity holds. |

### BAM QC And ID-SNP

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `BAM_QC:SAMBAMBA_PANEL_DEPTH` | `sambamba depth base` | Previously inside lowcov calculation | Internal only | Good tool split. |
| `BAM_QC:PANEL_DEPTH` | Python low-depth panel reducer | Old Perl `panel_depth.pl` | Internal only | Must preserve interval merge/order exactly. |
| `BAM_QC:LOWCOV_INTERSECT` | Bedtools `-loj` gene intersection | Previously inside lowcov flow | Internal only | Requires real BED file input, not a path string. |
| `BAM_QC:LOWCOV` | Python overlapping genes/lowcov finalization | Old Perl `overlapping_genes.pl`/LOWCOV | Publishes to `QC` | Highest-risk QC equivalence area. |
| `BAM_QC:BCFTOOLS_MPILEUP_CALL` | `bcftools mpileup` + call | First part of old `ALLELE_CALL` | Internal only | Good tool split. |
| `BAM_QC:BCFTOOLS_ANNOTATE_IDSNP` | `bcftools annotate` | Second part of old `ALLELE_CALL` | No QC publish rule | If `.final.vcf` should be visible in QC workflow, add matching QC publish rule. |
| `BAM_QC:BCFTOOLS_QUERY_IDSNP` | `bcftools query` | Third part of old `ALLELE_CALL` | Internal only | Feeds genotype JSON. |
| `BAM_QC:GENOTYPE2JSON` | Python genotype JSON | Old `genotype2json.py` | Publishes `*.genotypes.json` | Preserved. |
| `BAM_QC:SNP_CHECK` | Python ID-SNP concordance plus JSON combining | Old Perl `idsnp_controller.pl` plus shell commands | Publishes IDsnp JSONs | Current process expects `s${tumor}_c${normal}.json`; Python default preserves this if no `--out-prefix` changes it. |
| `BAM_QC:PAIRGEN_CDM` | Writes CDM load command | Same shell echo | Publishes `*.pairgen` | Preserved. |

### CNV Annotation

| Process path | Current command role | Old role | Published status | Notes |
|---|---|---|---|---|
| `CNV_ANNOTATE:COYOTE_SEGMENTS_BED` | Python creates raw segment BED | First part of old `coyote_segmentator.pl` | Internal only | Split required because Python no longer calls bedtools. |
| `CNV_ANNOTATE:COYOTE_SEGMENTS_INTERSECT` | Bedtools `-loj` with CNA genes | Subprocess inside old script | Internal only | Correct direction. |
| `CNV_ANNOTATE:COYOTE_SEGMENTS` | Python consumes intersect BED and writes final CN segment BEDs | Second part of old script | Publishes to `cnv` | Should preserve final `.cn-segments.*.bed` if script logic matches. |
| `CNV_ANNOTATE:GENE_INTERSECT` | Bedtools `-loj` before `cnvJSON.py` | Externalizes tool call | Internal only | Correct if input BED path is real and ordered. |
| `CNV_ANNOTATE:COYOTE_SEGMENTS_JSON` | Python JSON creation | Existing/new Python `cnvJSON.py` | Publishes panelmatched JSON | `cnvJSON.py` lacks explicit output args and derives names from `--id`; risk for config-driven naming. |
| `CNV_ANNOTATE:MERGE_SEGMENTS` | `cat` segment BEDs | Same simple merge | Publishes to `cnv` | OK, but order depends on channel list order. |
| `CNV_ANNOTATE:MERGE_JSON` | `jq` merge | Same merge | Publishes to `cnv` | OK for paired; tumor-only uses `cat`. |
| `CNV_ANNOTATE:SVDB_ANNOTATE_ARTEFACTS` | Shell helper wrapping multiple SVDB annotate DB passes | Old complex SVDB process | Publishes to `svvcf/merged` | Here-doc issue was fixed by moving DB logic out. Ensure helper emits exactly one `*.cnv.artefacts.vcf`. |

## Exact Result Equivalence Verdict

From static review, the branch is closer to a clean split but is not yet safe to declare exactly result-equivalent to `main`.

The major blockers are:

1. Several published intermediate files changed, especially GATK/CNVkit.
2. Some old published files are now only indirectly published or no longer published.
3. A few generic modules still have output-pattern risks (`FILTER_TNSCOPE`, `BEDTOOLS_INTERSECT`, `COYOTE_SEGMENTS_JSON`).
4. Python replacements need script-by-script parity confirmation, especially `find_contaminant.py`, `idsnp_controller.py`, `panel_depth.py`, `overlapping_genes.py`, `coyote_segmentator.py`, `mergeGATK*.py`, and `aggregate_vcf.py`.
5. Some config selectors are namespace-specific; a process may be correctly published under `ID_SNP:*` but not under `QC:*`.

## Removed And Added Process Names

Notable removed monolithic or misplaced process names:

- `modules/local/GATK/main.nf:FILTER_MERGE_GATK`
- `modules/local/GATK/main.nf:GATKCOV_CALL`
- `modules/local/GATK/main.nf:GATKCOV_COUNT`
- `modules/local/filters/main.nf:ANNOTATE_VEP`
- `modules/local/filters/main.nf:VCFANNO`
- `modules/local/filters/main.nf:CNV_BACKBONE_FILTER`
- `modules/local/filters/main.nf:FILTER_FOR_CNV`
- `modules/local/idSnp/main.nf:ALLELE_CALL`
- `modules/local/msisensor/main.nf:MSISENSOR`
- `modules/local/sentieon/main.nf:FILTER_TNSCOPE`
- `modules/local/sentieon/main.nf:SENTIEON_QC_TO_CDM`

Notable added split/tool-specific process names:

- `modules/local/bcftools/main.nf:ANNOTATE`
- `modules/local/bcftools/main.nf:BCFTOOLS_MPILEUP_CALL`
- `modules/local/bcftools/main.nf:BCFTOOLS_NORM`
- `modules/local/bcftools/main.nf:BCFTOOLS_QUERY_IDSNP`
- `modules/local/bedtools/main.nf:BEDTOOLS_INTERSECT`
- `modules/local/vcflib/main.nf:FILTER`
- `modules/local/vcflib/main.nf:GLXGT`
- `modules/local/vcftools/main.nf:VCF_CONCAT`
- `modules/local/vcftools/main.nf:VCF_SORT`
- `modules/local/tabix/main.nf:TABIX_BGZIPTABIX`
- `modules/local/vep/main.nf:ANNOTATE_VEP`
- `modules/local/vcfanno/main.nf:VCFANNO`
- `modules/local/pon/main.nf:CREATE_SNVPON`
- `modules/local/sambamba/main.nf:DEPTH`
- `modules/local/idSnp/main.nf:GENOTYPE2JSON`
- `modules/local/qc/main.nf:PANEL_DEPTH`
- `modules/local/qc/main.nf:CONTAMINATION_CDM`
- Multiple split GATK CNV processes

This direction is consistent with the tool-based module split, but channel and output equivalence still need runtime verification.

## Recommendation

Do not claim exact same results yet.

Recommended next audit steps before testing full production profiles:

1. Build a file-level expected-output matrix from `main` for one tumor-only and one tumor-normal case.
2. Compare current branch process publish rules against that matrix.
3. Run a minimal `-stub-run` only after syntax review, then run one small real case per profile.
4. Compare checksums or normalized content for final VCFs, CNV JSON, contamination, idSNP, CNVkit GENS, and Coyote upload files.

The migration is directionally coherent, but several output capture and result-equivalence risks remain.

---

# Command-Level Split Audit

This section answers the more specific question: when `main` had several commands inside one process, did the `containers` branch split those commands into separate processes while preserving the produced final files?

This still does not prove byte-identical runtime results, but it is a closer static check than channel-shape comparison.

## Summary By Area

| Area | Old monolithic process | Current split | Command sequence preserved? | Final file capture status |
| --- | --- | --- | --- | --- |
| FreeBayes SNV calling | `FREEBAYES` | `FREEBAYES -> FREEBAYES_FILTER_LOWCOV -> FREEBAYES_FILTER_LOWFRQ -> FREEBAYES_VCFGLXGT -> FILTER_FREEBAYES -> FREEBAYES_REMOVE_AD` | Mostly yes | Likely preserved, check final `freebayes_${bed}.vcf` naming |
| VarDict SNV calling | `VARDICT` plus filter in same/broader module path | `VARDICT -> FILTER_VARDICT` | Likely yes | Needs content parity check against Perl filter |
| TNScope SNV calling | `TNSCOPE` and `FILTER_TNSCOPE` in Sentieon module | `TNSCOPE -> FILTER_TNSCOPE` | Likely yes | Needs content parity check against Perl filter |
| Pindel | `PINDEL_CALL` doing call and conversion | `PINDEL_CALL -> PINDEL_TO_VCF -> FILTER_PINDEL` | Mostly yes | Likely preserved if final `*_pindel.vcf` matches old |
| CNVkit call/export | `CNVKIT_CALL` ran call, export VCF, and nexus-ogt | `CNVKIT_CALL -> CNVKIT_EXPORT_VCF -> CNVKIT_EXPORT_NEXUS_OGT` | Yes | VCF and `.cnvkit` files still produced, but publish ownership changed |
| CNVkit GENS | `CNVKIT_GENS` ran Perl and produced `.bed.gz`; `MERGE_GENS` handled tabix/merge | `CNVKIT_GENS -> TABIX_BGZIPTABIX -> MERGE_GENS` | Mostly yes | Same final `.bed.gz*`, `.gens`, `.gens_v4_somatic` intended; verify file names |
| GATK tumor CNV | `GATKCOV_COUNT`, `GATKCOV_CALL`, `GATK2VCF`, `MERGE_GATK_TUMOR` | many single-command GATK processes plus Python merge | Mostly yes | Final `*_gatk_tumor_merged.vcf` preserved; intermediate published set changed |
| GATK germline CNV | `GATK_COUNT_GERMLINE -> GATK_CALL_PLOIDY -> GATK_CALL_GERMLINE_CNV -> POSTPROCESS -> FILTER_MERGE_GATK` | same first three, plus `GATK_EXTRACT_GERMLINE_CNV_SHARDS -> POSTPROCESS -> FILTER_GATK -> MERGE_GATK` | Mostly yes | Final `*.gatk.filtered.merged.vcf` preserved; shard extraction split needs runtime check |
| BAM QC low coverage | `LOWCOV` ran `panel_depth.pl` and `overlapping_genes.pl` | `SAMBAMBA_PANEL_DEPTH -> PANEL_DEPTH -> BEDTOOLS_INTERSECT -> LOWCOV` | Conceptually yes, implementation changed more than split | High risk: output equivalence depends on Python replacements and explicit bedtools step |
| idSNP | `ALLELE_CALL` ran four commands | `BCFTOOLS_MPILEUP_CALL -> BCFTOOLS_ANNOTATE_IDSNP -> BCFTOOLS_QUERY_IDSNP -> GENOTYPE2JSON` | Yes | Final `.final.vcf` and `.genotypes.json` intended; publish rules look preserved for JSON only |
| Contamination | `CONTAMINATION` ran `find_contaminant.pl`, echo, paste | `CONTAMINATION -> CONTAMINATION_CDM` | Yes structurally | Must verify Python produces same `.txt`, `.png`, `.value`; `.contaminationpy` moved to second process |

## FreeBayes

Old `main` process `FREEBAYES` ran these commands:

1. `freebayes ... > freebayes_${bed}.vcf.raw`
2. `vcffilter ... | vcffilter ... | vcfglxgt > freebayes_${bed}.filt1.vcf`
3. `filter_freebayes_somatic.pl` or `filter_freebayes_unpaired.pl > freebayes_filtered_${bed}.vcf`
4. `bcftools annotate -x FORMAT/AD ... -o freebayes_${bed}.vcf`

Current split:

1. `FREEBAYES` runs `freebayes ... > freebayes_${bed}.vcf.raw`
2. `FREEBAYES_FILTER_LOWCOV` runs first `vcffilter`
3. `FREEBAYES_FILTER_LOWFRQ` runs second `vcffilter`
4. `FREEBAYES_VCFGLXGT` runs `vcfglxgt`
5. `FILTER_FREEBAYES` runs Python replacement for the Perl somatic/unpaired filter
6. `FREEBAYES_REMOVE_AD` runs `bcftools annotate -x FORMAT/AD`

Assessment: command split is structurally correct. The final expected file should still be `freebayes_${bed}.vcf`. The main risk is whether `FILTER_FREEBAYES` Python is exactly equivalent to both Perl branches and whether `ext.prefix`/`ext.suffix` reconstructs the same final basename for every bed part.

## CNVkit

Old `CNVKIT_CALL` ran three independent outputs in one process:

1. `cnvkit.py call`
2. `cnvkit.py export vcf`
3. `cnvkit.py export nexus-ogt`

Current split:

1. `CNVKIT_CALL` only runs `cnvkit.py call`
2. `CNVKIT_EXPORT_VCF` runs `cnvkit.py export vcf`
3. `CNVKIT_EXPORT_NEXUS_OGT` runs `cnvkit.py export nexus-ogt`

Assessment: this is the intended split. Final files are still represented:

- `*.${part}.call*.cns`
- `*.${part}.cnvkit.vcf`
- `*.${part}_logr_ballele.cnvkit`

Publish behavior changed: old `CNVKIT_CALL` published call CNS and VCF together; now `CNVKIT_CALL`, `CNVKIT_EXPORT_VCF`, and `CNVKIT_EXPORT_NEXUS_OGT` publish separately. That is fine if the desired final published file set includes all three. If `.cnvkit` files were not previously published intentionally, this is an added output.

## CNVkit GENS

Old `CNVKIT_GENS`:

1. `generate_gens_data_from_cnvkit.pl`
2. `mv ${meta.id}.baf.bed.gz ${prefix}.${part}.baf.bed.gz`
3. `mv ${meta.id}.cov.bed.gz ${prefix}.${part}.cov.bed.gz`

Current split:

1. `CNVKIT_GENS` runs `generate_gens_data_from_cnvkit.py` and produces uncompressed `.baf.bed` and `.cov.bed`
2. `CNVKIT_BAF_BGZIPTABIX` compresses/indexes BAF bed
3. `CNVKIT_COV_BGZIPTABIX` compresses/indexes coverage bed
4. `MERGE_GENS` consumes `.bed.gz` files and creates merged `.bed.gz`, `.tbi`, `.gens`, and `.gens_v4_somatic`

Assessment: command responsibility is correctly split by tool, but this is not a simple mechanical split because the Python script now emits uncompressed beds and compression moved out. Final files appear intended to match:

- `${prefix}.${part}.baf.bed.gz`
- `${prefix}.${part}.cov.bed.gz`
- `${prefix}.merged.sorted.baf.bed.gz`
- `${prefix}.merged.sorted.cov.bed.gz`
- `${prefix}.gens`
- `${prefix}.gens_v4_somatic`

Risk: `MERGE_GENS` publishes `*.bed.gz*` from its work directory, not directly from `CNVKIT_GENS`. This can preserve output capture only if the compressed files are staged into `MERGE_GENS` with the same names.

## GATK CNV

Old tumor command grouping:

- `GATKCOV_COUNT` ran `CollectReadCounts`, `DenoiseReadCounts`, and `PlotDenoisedCopyRatios`
- `GATKCOV_CALL` ran `ModelSegments`, `CallCopyRatioSegments`, and `PlotModeledSegments`
- `GATK2VCF` ran `gatk_to_vcf.py`
- `MERGE_GATK_TUMOR` ran `mergeGATK_tumor.pl`

Current split:

- `GATK_COLLECT_READ_COUNTS`
- `GATK_DENOISE_READ_COUNTS`
- `GATK_PLOT_DENOISED_COPY_RATIOS`
- `GATK_MODEL_SEGMENTS`
- `GATK_CALL_COPY_RATIO_SEGMENTS`
- `GATK_PLOT_MODELED_SEGMENTS`
- `GATK2VCF`
- `MERGE_GATK_TUMOR`

Assessment: command sequence is preserved and split well. Final tumor file `*_gatk_tumor_merged.vcf` is still produced. The old Perl merge is now Python, so exact content depends on Python parity.

Old germline command grouping:

- `GATK_COUNT_GERMLINE`
- `GATK_CALL_PLOIDY`
- `GATK_CALL_GERMLINE_CNV`
- `POSTPROCESS` extracted tar files and ran `PostprocessGermlineCNVCalls`
- `FILTER_MERGE_GATK` ran `filter_gatk.pl` and `mergeGATK.pl`

Current split:

- same first three commands
- `GATK_EXTRACT_GERMLINE_CNV_SHARDS` extracts tar/ploidy directories
- `POSTPROCESS` only runs `PostprocessGermlineCNVCalls`
- `FILTER_GATK` runs Python replacement for `filter_gatk.pl`
- `MERGE_GATK` runs Python replacement for `mergeGATK.pl`

Assessment: command sequence is preserved, and this is a clean split. Risk is lower than before, but exact results still require parity checks for the Python replacements and tar extraction layout.

## BAM QC Low Coverage

Old `LOWCOV` process:

1. `panel_depth.pl $bam $args > lowcov.bed`
2. `overlapping_genes.pl lowcov.bed $args2 > ${prefix}.lowcov.bed`

Current split:

1. `SAMBAMBA_PANEL_DEPTH` runs `sambamba depth`
2. `PANEL_DEPTH` runs `panel_depth.py`
3. `LOWCOV_INTERSECT` runs `bedtools intersect -loj`
4. `LOWCOV` runs `overlapping_genes.py`

Assessment: this is more than just splitting existing shell commands; the algorithmic boundary changed because `bedtools intersect` is now an explicit process and both Perl scripts were replaced. The final file name `${prefix}.lowcov.bed` is still intended. This area needs golden-output comparison.

## idSNP

Old `ALLELE_CALL` process:

1. `bcftools mpileup | bcftools call > ${prefix}.raw.vcf`
2. `bcftools annotate ... -o ${prefix}.final.vcf`
3. `bcftools query ... > ${prefix}.genotypes`
4. `genotype2json.py ... ${prefix}.genotypes.json`

Current split:

1. `BCFTOOLS_MPILEUP_CALL`
2. `BCFTOOLS_ANNOTATE_IDSNP`
3. `BCFTOOLS_QUERY_IDSNP`
4. `GENOTYPE2JSON`

Assessment: command split is correct. Final `.final.vcf` and `.genotypes.json` are still produced. Publish behavior should be reviewed: current config publishes `.genotypes.json`; old config also only kept `.genotypes.json` from this group while discarding raw/final VCF except as internal process inputs.

`SNP_CHECK` is not split by command; it replaces `idsnp_controller.pl` with `idsnp_controller.py` and keeps the same downstream `cp`, `echo`, and `combinejsons.py` flow. Final `.T.idsnp`, `.N.idsnp`, and `.pairgen` names are preserved.

## Contamination

Old `CONTAMINATION` process:

1. `find_contaminant.pl ... > ${sample}.value`
2. `echo ... > ${sample}.1`
3. `paste ... > ${sample}.contaminationpy`

Current split:

1. `CONTAMINATION` runs `find_contaminant.py ... --out ${sample}.value`
2. `CONTAMINATION_CDM` runs `echo` and `paste` to produce `.contaminationpy`

Assessment: command split is correct in shape. The final `.contaminationpy` moved to a separate process and should still be published by `CONTAMINATION_CDM`. The important parity question is whether `find_contaminant.py` always creates the same `.txt` and `.png` side products as `find_contaminant.pl`, because `CONTAMINATION` still declares those as outputs.

## Revised Conclusion

The split itself is mostly aligned with the intended rule: separate tools/commands are now separate processes.

The remaining question is not channel names. It is final-file parity. From static review:

- FreeBayes split: mostly correct, final filename likely preserved.
- CNVkit call/export split: correct, with changed publish ownership.
- CNVkit GENS split: conceptually correct, but final `.bed.gz*` capture depends on downstream staging through `MERGE_GENS`.
- GATK split: mostly correct and cleaner than before.
- idSNP split: command order preserved.
- Contamination split: command order preserved, but Python side products must match Perl.
- Lowcov split: highest risk because it is a rewritten flow, not only a split.

So the next report/testing task should be a final-file manifest comparison, not another channel-shape audit.
