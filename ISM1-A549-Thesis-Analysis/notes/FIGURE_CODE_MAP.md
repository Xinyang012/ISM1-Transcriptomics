# Preliminary figure-to-script map (requires validation)

- **Figure 2 / single-cell panels:** `03_scRNAseq_GSE292515/01_import_qc_umap_original.R` contains the Ism1 FeaturePlot for young/old Mock/IAV and the QC processing. Confirm filenames with the exact figure assembly.
- **Figure 2 / public bulk panels:** `02_Public_Bulk_GSE271812/` is incomplete. Original upstream scripts are not yet recovered.
- **Figure 3 / laboratory bulk RNA-seq:** `01_Lab_Bulk_RNAseq/analysis_draft_from_original.R` contains laboratory differential expression and GSEA functions. Match every panel with the exact output file.
- **Supplementary figures:** exploration scripts contain additional candidate plots, but mapping is not verified.

- **Public bulk GSE271812:** `02_Public_Bulk_GSE271812/GSE271812_original_downstream_analysis_REVIEW.R` includes Ism1–gene correlation, GSEA/GSVA and plotting explorations; requires missing upstream R objects; final Figure 2 panels must be matched before release.
