# Validation and thesis cross-checks

- [ ] Confirm sample metadata mapping (3 h/12 h; NOR/SEN; +/- rISM1) and sequencing replicates against lab records.
- [ ] Confirm final DESeq2 contrasts, adjusted p-value and |log2FC| cutoffs for each reported heatmap. Do **not** assume every plot used the same cutoff.
- [ ] Recover the complete public bulk GSE271812 analysis script, not only correlation fragments.
- [ ] Verify 10 mouse samples' age/treatment/replicate labels for GSE292515.
- [ ] Confirm exact GSE292515 QC counts before and after filtering; check excluded-cell patterns per sample.
- [ ] Confirm how Seurat v5 layers were joined, Harmony run parameters and identity assignment in the exported plots.
- [ ] Confirm the published Figure 2 visualisation uses the documented Ism1 FeaturePlot object and colour scale.
- [ ] Verify no trial-only IFN module score or exploratory pseudobulk result is misrepresented as a thesis figure.
- [ ] Run all retained final scripts locally and save `sessionInfo()` output; R execution not available in this build environment.
- [ ] Supervisor approval for release of unpublished lab analysis scripts and any metadata.
- [ ] If repository remains private, ensure examiners can access it before citing the repository link.

- [ ] Recover upstream GSE271812 object creation (`vst.matrix.iav271812`, `metadata.timecourse271812`, `count.matrix.primary271812`, `gsva.plot.data.271812`) and test from fresh R session.
- [ ] Verify GSE271812 Figure 2 panel mapping, statistical cutoff, group labels, and complete input metadata.
