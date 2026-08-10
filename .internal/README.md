# Internal Published Output Parity Checks

These files are for validating the containers refactor against `main`/`master`.
They compare only published output directories, not Nextflow work directories or
unpublished intermediates.

## Files

- `published_output_test_cases.json`: paired/unpaired profile coverage and command metadata.
- `render_test_commands.py`: renders baseline and candidate Nextflow command templates.
- `build_published_manifest.py`: records relative path, size, and SHA-256 for every published file.
- `compare_published_manifests.py`: compares baseline and candidate manifests.

## Workflow

1. Render the commands:

   ```bash
   .internal/render_test_commands.py
   ```

2. Run each baseline command on `main`/`master` with `BASELINE_OUTDIR` set.
3. Run the matching candidate command on this branch with `CANDIDATE_OUTDIR` set.
4. Build manifests for the published subdirectory from each run:

   ```bash
   .internal/build_published_manifest.py \
     --case-id spp_solid_paired \
     --published-dir "$BASELINE_OUTDIR/solid_hg38" \
     --out .internal/results/spp_solid_paired.baseline.json

   .internal/build_published_manifest.py \
     --case-id spp_solid_paired \
     --published-dir "$CANDIDATE_OUTDIR/solid_hg38" \
     --out .internal/results/spp_solid_paired.candidate.json
   ```

5. Compare the manifests:

   ```bash
   .internal/compare_published_manifests.py \
     --baseline .internal/results/spp_solid_paired.baseline.json \
     --candidate .internal/results/spp_solid_paired.candidate.json \
     --json-out .internal/results/spp_solid_paired.compare.json \
     --md-out .internal/results/spp_solid_paired.compare.md
   ```

The comparator exits with status `1` when any published file is missing, extra,
or content-different.
