<!-- Pull Request Template for SomaticPanelPipeline -->

# Summary
Clearly describe the purpose of this PR and what changed. 
If it fixes an issue, reference it (e.g., Fixes #123). 
For reviewers: include how to test this change.  

---

## Type of change
- [ ] Bug fix  
- [ ] Patch / hotfix  
- [ ] New feature / enhancement  
- [ ] Refactor / cleanup  
- [ ] Documentation update  


---

## Checklist (author)  

- [ ] **CHANGELOG** updated with a clear entry   
- [ ] **Version bumped** in `configs/nextflow.hopper.config` (if user-facing change, skip if it is an infra only)  
- [ ] Pipeline runs successfully on stubs (`-profile test` or equivalent)  
- [ ] Outputs validated (VCFs, metrics, reports, MultiQC if applicable)  
- [ ] Output from affected processes inspected (`.command.sh`, `.command.log`, `.command.err`, result files in `work/`)  
- [ ] New data formats tested for compatibility with **Coyote3**  
- [ ] New data formats tested for compatibility with **GENS**  
- [ ] Documentation updated (README, params, examples, usage)  
- [ ] Containers / tools updated are version-pinned and reproducible  
- [ ] At least one reviewer has tested and approved the code  

---

## Breaking changes
- [ ] No breaking changes  
- [ ] Params/outputs/configs changed (list details in **Summary**)  
- [ ] Output filenames or structure changed (document migration path)  

---

## Testing performed
Describe what was tested, datasets/profiles used, and results. Attach logs or links if helpful.

Performed by:  
- [ ] Viktor  
- [ ] Ram  
- [ ] Sailedra  

(Add if missing)

---

## Notes for reviewers
Mention risk areas, edge cases, or anything needing special attention. 

---

## Review performed by
- [ ] Viktor  
- [ ] Ram  
- [ ] Sailedra  

(Add if missing)

---

## Release notes (draft)
Short, user-facing bullet points for the next release tag.
