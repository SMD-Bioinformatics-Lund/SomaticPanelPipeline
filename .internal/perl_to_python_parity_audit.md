# Perl To Python Parity Audit

Date: 2026-08-10

Scope: static review of legacy Perl scripts from `main:bin/` against Python replacements in the current `containers` working tree. I did not execute any scripts, run containers, or compare real output files.

The goal of this audit is not to require identical function names or implementation style. The requirement is output-contract parity: same intended files, same record logic, same order where order was meaningful, same filtering decisions, and explicit file outputs instead of relying on stdout where the refactor is improving script behavior.

## Summary

The Python migration is mostly in the right direction: most old stdout-based Perl scripts now have `argparse` and explicit `--out` output paths. Several scripts use shared VCF helpers, which is better than copy-pasted parsing.

However, exact output equivalence is not proven. The highest-risk replacements are scripts that either:

- produce multiple side-effect files,
- previously called external tools internally,
- affect downstream database/cron files,
- or are sensitive to row order/formatting.

Highest-risk scripts:

- `find_contaminant.py`
- `idsnp_controller.py`
- `panel_depth.py`
- `overlapping_genes.py`
- `coyote_segmentator.py`
- `cnvJSON.py`
- `aggregate_vcf.py`
- `mergeGATK.py`
- `mergeGATK_tumor.py`
- `generate_gens_data_from_cnvkit.py`

## Scripts Without `argparse`

These current Python files do not have an `argparse` CLI. They are acceptable if they are library/helper modules only and are not called directly by Nextflow.

| File | Status | Action |
|---|---|---|
| `bin/cmdvcf.py` | Helper/library style | OK if never called as a process command. |
| `bin/vcf2.py` | Replacement parser/helper for legacy `vcf2.pm` | OK as a library. |
| `bin/vcf_pipeline_utils.py` | Shared helper functions | OK as a library. |

If any of these are called directly from a process, they should get a CLI parser or be moved clearly under a library/helper folder.

## Scripts With CLI But Not Explicit `--out`

These scripts have a parser but still rely on positional outputs, output prefixes, or derived output names. This has now been tightened for the active process-called scripts: the Python CLIs keep their previous positional/prefix behavior for compatibility, but the Nextflow process commands pass explicit output paths where those scripts are used.

| File | Current output style | Risk | Recommendation |
|---|---|---|---|
| `bin/cnvJSON.py` | Supports `--out`, `--panelmatched-out`, and compatibility `--out-prefix`/`--id` derived names | Low. Active `COYOTE_SEGMENTS_JSON` now passes both output files explicitly | Fixed. |
| `bin/coyote_segmentator.py` | Supports explicit `--segments-out` for interval extraction and `--panel-out`/`--raw-out` for final BEDs | Low. Active processes now pass explicit BED outputs | Fixed. |
| `bin/filter_manta.py` | Supports explicit `--filtered-out` and `--bnd-out`, with compatibility `--out-prefix` | Low. Active `FILTER_MANTA*` process command passes both VCF paths explicitly | Fixed. |
| `bin/generate_gens_data_from_cnvkit.py` | Supports explicit `--baf-out` and `--cov-out`, with compatibility `--out-prefix` | Low. Active `CNVKIT_GENS` process now passes both BED paths explicitly | Fixed. |
| `bin/idsnp_controller.py` | Supports explicit `--csv-out` and `--json-out`, with compatibility `--out-prefix` | Low/medium. Active `SNP_CHECK` now passes the exact legacy `s<TUMOR>_c<NORMAL>` CSV/JSON names explicitly | Fixed, but sample/control ordering should still be golden-tested. |
| `bin/genotype2json.py` | Supports explicit `--input` and `--out`, with positional compatibility | Low. Active `GENOTYPE2JSON` now uses named args | Fixed. |
| `bin/combinejsons.py` | Supports explicit `--info-json`, `--genotype-json`, `--partner-run-json-file`, and `--out`, with positional and old underscore-option compatibility | Low. Active `SNP_CHECK` now uses named args | Fixed. |
| `bin/compare_vcf.py` | Positional VCFs, output prefix | Not part of old Perl migration and likely diagnostic | OK if not used in production workflow. |

## Scripts Still Writing To Stdout

Most current scripts write main data to explicit files. The remaining `print(..., file=out_fh)` calls are data writes to opened output files and are acceptable.

Cleanup completed:

- `combinejsons.py` no longer prints a status line to stdout.
- `idsnp_controller.py` no longer prints a status line to stdout.
- `cnvJSON.py` now reports file-not-found errors to stderr.

Remaining intentional stdout cases:

- `compare_vcf.py` intentionally prints a human-readable comparison summary and appears diagnostic rather than production workflow code.
- `vcf_pipeline_utils.open_output()` and a few VCF filters still support `--out -` for manual CLI compatibility. Active Nextflow process commands pass real output files, so production workflow output is not stdout-driven.
- `post_annotation_filtering.py` still supports `--out -`, but its process command passes `--out <file>`.

The old Perl pattern of `script.pl ... > output.file` has been replaced for the active Python migrations reviewed here. Any remaining shell redirection in module scripts is for non-Python tools, such as commands that naturally write VCF to stdout.

## Perl To Python Mapping

| Legacy Perl | Current Python | Old output behavior | Current output behavior | Parity assessment |
|---|---|---|---|---|
| `aggregate_vcf.pl` | `aggregate_vcf.py` | Wrote aggregate VCF, likely stdout or process-managed output | Uses `--vcf`, `--sample-order`, `--out` | High risk. Aggregation order, caller priority, sample column order, FILTER/header order, and duplicate handling must be golden-tested. |
| `cnvkit2HRD.pl` | `cnvkit2HRD.py` | ARGV/stdout-style HRD table conversion | `--input`, `--sample-id`, `--ploidy`, `--out` | Likely preserved. Check exact column order and numeric formatting. |
| `coyote_segmentator.pl` | `coyote_segmentator.py` | One script parsed VCF, ran bedtools intersect internally, wrote segment BED outputs | Split into raw BED creation, `BEDTOOLS_INTERSECT`, and Python finalization. Process commands pass `--segments-out`, `--panel-out`, and `--raw-out` explicitly | Structurally improved. Exact parity depends on external bedtools command matching old `-loj` call and preserving sorted input order. |
| `create_snv_pon.pl` | `create_snv_pon.py` | ARGV/stdout PON table | `--vcf-mask`, `--out` | Likely preserved. Needs check for glob ordering and table sort order. |
| `filter_cnvkit.pl` | `filter_cnvkit.py` | ARGV/stdout filter | `--tsv`, `--coverage`, `--distance`, `--out` | Likely preserved. Check float/int comparisons and header retention. |
| `filter_freebayes_somatic.pl` | `filter_freebayes_somatic.py` | Perl VCF filter to stdout | `--vcf`, `--tumor`, `--normal`, `--out`, threshold args | Medium risk. Core thresholds appear preserved; verify INFO/FORMAT field parsing, filter order, and final trailing fields. |
| `filter_freebayes_unpaired.pl` | `filter_freebayes_unpaired.py` | Perl VCF filter to stdout | `--vcf`, `--out`, threshold args | Medium risk. Header quirks and WARN/FAIL behavior should be compared. |
| `filter_gatk.pl` | `filter_gatk.py` | Perl GATK VCF filter to stdout | `--vcf`, `--out` | Medium risk. Header injection and record rewrite should be compared. |
| `filter_manta.pl` | `filter_manta.py` | Created two files: `_manta_filtered.vcf` and `_manta_bnd_filtered.vcf` | Uses `--vcf`, `--sample-id/--id`, `--af`, and explicit `--filtered-out`/`--bnd-out` in active process commands; `--out-prefix` remains for compatibility | Likely preserved. Current Nextflow no longer relies on default-name drift. |
| `filter_pindel_somatic.pl` | `filter_pindel_somatic.py` | Perl VCF filter to stdout | `--vcf`, `--out` | Likely preserved. Low complexity. |
| `filter_tnscope_somatic.pl` | `filter_tnscope_somatic.py` | Perl paired TNscope VCF filter to stdout | `--vcf`, `--tumor`, `--normal`, `--out`, threshold args | Medium risk. Filter decision parity needs golden examples. |
| `filter_tnscope_unpaired.pl` | `filter_tnscope_unpaired.py` | Perl unpaired TNscope VCF filter to stdout | `--vcf`, `--out`, threshold args | Medium risk. |
| `filter_vardict_somatic.pl` | `filter_vardict_somatic.py` | Perl paired VarDict VCF filter to stdout | `--vcf`, `--tumor`, `--normal`, `--out`, threshold args | Medium risk. Core thresholds appear represented. |
| `filter_vardict_unpaired.pl` | `filter_vardict_unpaired.py` | Perl unpaired VarDict VCF filter to stdout | `--vcf`, `--out`, threshold args | Medium risk. |
| `filter_with_ffpe_pon.pl` | `filter_with_ffpe_pon.py` | Perl PON/FFPE filter to stdout | `--vcf`, `--tumor-id`, `--pons`, `--out` | Medium risk. Must compare caller-specific PON fields and appended FILTER values. |
| `filter_with_pon.pl` | `filter_with_pon.py` | Perl PON filter to stdout | `--vcf`, `--tumor-id`, `--pons`, `--out` | Medium risk. |
| `find_contaminant.pl` | `find_contaminant.py` | Printed contamination value to stdout and wrote side files `${id}.dist.txt`, `${id}.genotypes.txt`, `${id}.png` | Writes value via `--out`; still creates side files from sample id/paired id | High risk. Logic is largely mirrored, but PNG rendering differs from GD::Graph and no-bin side-file behavior is not identical: old touched PNG only; current also touches dist file. |
| `fix_vep_gnomad.pl` | `fix_vep_gnomad.py` | Perl VEP CSQ field normalization to stdout | `--vcf`, `--out` | Likely preserved, but CSQ allele ordering needs targeted examples. |
| `gatksegments.pl` | `gatksegments.py` | ARGV/stdout segment table | `--cr`, `--baf` accepted/ignored, `--purity`, `--out` | Likely preserved if BAF was also ignored or unused before. If Perl used BAF, Python is not equivalent. |
| `generate_gens_data_from_cnvkit.pl` | `generate_gens_data_from_cnvkit.py` | Generated BAF/COV BED and compressed/indexed outputs in one flow | Active process passes explicit `--baf-out` and `--cov-out`; compression moved to tabix processes | Content may be preserved, but output lifecycle differs. Need compare pre-compression BED content and final `.gz/.tbi` names. |
| `germline_for_cnvkit.pl` | `germline_for_cnvkit.py` | Perl VCF germline conversion to stdout | `--vcf`, `--out`, threshold args | Medium risk. It also contains mark-germline helper code, so ensure the called CLI path is the intended one. |
| `idsnp_controller.pl` | `idsnp_controller.py` | Created `s<SAMPLE>_c<CONTROL>.csv` and `.json`; printed usage/status | Active process passes explicit `--csv-out s<TUMOR>_c<NORMAL>.csv` and `--json-out s<TUMOR>_c<NORMAL>.json`; compatibility default remains | Medium/high risk. The filename is no longer implicit, but sample/control ordering and concordance metrics still need golden comparison. |
| `mark_germlines.pl` | `mark_germlines.py` | Perl VCF germline marking to stdout | `--vcf`, `--tumor-id`, optional `--normal-id`, `--assay`, `--out` | Medium risk. Check paired/unpaired output filename and FILTER/INFO injection. |
| `mergeGATK.pl` | `mergeGATK.py` | Perl GATK normal merge to stdout/file | `--vcf`, `--out` | High risk. Merging adjacent/overlapping events is order-sensitive. |
| `mergeGATK_tumor.pl` | `mergeGATK_tumor.py` | Perl tumor GATK merge | `--vcf`, `--out` | High risk. Earlier copied-run error was script invocation/shebang related, not necessarily logic, but content parity still needs comparison. |
| `merge_melt.pl` | `merge_melt.py` | Perl MELT merge, including empty-header behavior | `--vcf-header`, `--sample-id`, `--out` | Medium risk. Empty input/header case should be golden-tested. |
| `overlapping_genes.pl` | `overlapping_genes.py` | Perl lowcov gene overlap summarizer | Uses parser with input/output args | High risk. This is a QC-reporting output where exact row grouping/order matters. |
| `panel_depth.pl` | `panel_depth.py` | Perl panel-depth reduction, likely stdout/file via process | `--depth`, `--cutoff`, `--out` | High risk. Numeric cutoff and interval merge behavior must match exactly. |
| `provider.pl` | `provider.py` | Perl provider fingerprint generator, also ran pileup-related external work | `--bam`, `--out`, `--pileup`, `--bed`, `--bedxy`, flags | Medium/high risk if used. It still has external-tool design considerations unless pileup is always precomputed. |
| `qc_sentieon.pl` | `qc_sentieon.py` | Perl JSON/QC command output | `--sample-id`, `--type`, `--out` | Medium risk. Check JSON keys, numeric formatting, and missing-value handling. |

## External Tool Calls Inside Python

The refactor goal is that Python scripts should not call separate bioinformatics tools like bedtools/samtools/d4tools when those can be represented as Nextflow processes.

Current static status:

| Script | External command use | Status |
|---|---|---|
| `coyote_segmentator.py` | Previously had bedtools-style functionality; current CLI supports `--segments-out` and `--intersect-bed` | Good direction if the subprocess call is fully removed from active path. The workflow now has `COYOTE_SEGMENTS_INTERSECT`. |
| `coyote_d4_cov.py` | Uses `subprocess` for `d4tools`/coverage work | Exception accepted by user because D4 uses a dedicated container. |
| `provider.py` | Has logic around BAM/pileup generation and may need external samtools workflow separation depending active use | Should be reviewed if provider is used by any active workflow. |
| `idsnp_controller.py` | No external tool subprocess; consumes prepared VCFs | Good. |
| `panel_depth.py`, `overlapping_genes.py` | No external tool subprocess in the intended split; consume sambamba/bedtools outputs | Good if workflows always supply precomputed inputs. |

## Output Contract Details To Golden-Test

The following fixture comparisons should be done before declaring parity:

1. FreeBayes paired and tumor-only VCF:
   - old `FREEBAYES` final filtered VCF after LowCov, LowFrq, `vcfglxgt`, Perl filter, and AD removal
   - current split output after `FREEBAYES_REMOVE_AD`

2. VarDict paired and tumor-only VCF:
   - old Perl filter output
   - current `filter_vardict_*.py` output

3. TNscope paired and tumor-only VCF:
   - old Perl filter output
   - current `filter_tnscope_*.py` output

4. PON and FFPE PON:
   - verify identical FILTER fields, header additions, and caller-specific annotation behavior.

5. Contamination:
   - tumor-only and tumor-normal cases
   - exact `.value` text
   - expected presence of `.dist.txt`, `.genotypes.txt`, `.png`
   - normal-sample side-file prefix when `--normal` is used

6. ID-SNP:
   - exact `s<TUMOR>_c<NORMAL>.csv`
   - exact `s<TUMOR>_c<NORMAL>.json`
   - final `<sample>.T.idsnp` and `<normal>.N.idsnp`
   - behavior when sample/control order is reversed upstream

7. CNVkit GENS:
   - exact `.baf.bed` and `.cov.bed` content before compression
   - exact `.bed.gz`, `.bed.gz.tbi`, `.gens`, and `.gens_v4_somatic` names

8. CNV annotation:
   - `coyote_segmentator.py` raw and panel BEDs
   - `cnvJSON.py` JSON fields and `nprobes` numeric handling
   - merge order for paired JSONs and segment BEDs

9. GATK merge/filter:
   - `filter_gatk.py`
   - `mergeGATK.py`
   - `mergeGATK_tumor.py`
   - `gatk_to_vcf.py`

## Recommended Python CLI Standard

For new and migrated scripts:

1. Every process-called Python script should use `argparse`.
2. Single-output scripts should require `--out`.
3. Multi-output scripts should either:
   - require explicit output paths for each file, or
   - use `--out-prefix` and document the exact derived filenames in the module.
4. Scripts should write data to files, not stdout.
5. Status messages should go to stderr or be removed.
6. Python should not call other bioinformatics tools unless that tool cannot reasonably be split into a Nextflow process. The accepted current exception is D4-related code.

## Verdict

The Python conversions are directionally correct but not yet proven exact. Most scripts have the right CLI shape now. The remaining work is not broad rewrites; it is targeted parity hardening:

- make all output paths explicit where names matter,
- remove or document derived-name side effects,
- ensure subprocess tool calls are represented as Nextflow processes,
- and run golden-output comparisons when testing is allowed.
