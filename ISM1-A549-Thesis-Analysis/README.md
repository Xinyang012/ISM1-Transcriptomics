# ISM1-A549 Thesis Analysis — working code archive

This repository is a **pre-publication, unverified collection of the user's existing R analyses**, organised for MSc Biotechnology thesis recordkeeping. **Not a polished or fully reproducible release yet.** The laboratory data, GEO inputs, R objects and package environment are not included.

## Contents

| Location | Analysis | Status |
|---|---|---|
| `01_Lab_Bulk_RNAseq/` | A549 3 h / 12 h RNA-seq, DESeq2, Hallmark GSEA and plotting | Original analysis script, needs input tables and validation |
| `02_Public_Bulk_GSE271812/` | GSE271812 public bulk RNA-seq, Spearman/GO-KEGG GSEA/GSVA | **Recovered substantial original downstream analysis script**; missing upstream data/object initialization; not yet standalone |
| `03_scRNAseq_GSE292515/01_import_qc_umap_original.R` | Download/import, QC, Seurat, UMAP and Ism1 expression display | Original draft, repeated commands and environment-specific paths may remain |
| `03_scRNAseq_GSE292515/02_downstream_exploratory_original.R` | Additional clustering, pseudobulk, ISG, and exploratory analyses | Original draft; some analyses may not appear in final figures |

## Important limitations

- Original code was preserved as closely as practical; the user's local absolute path was replaced with `PATH_TO_*` placeholder strings. Replace them on a local machine before execution.
- This archive has **not** been run against the underlying count matrices or Seurat objects; do not claim that all figures can be regenerated yet.
- No original sequencing data or internal working files are included. Check data-sharing permission with the supervisor.
- Not every exploratory analysis in `02_downstream_exploratory_original.R` was used in the thesis. Retain or move these after checking the figure provenance.
- Original scripts contain repeated code or interactive package installation commands. Leave them as a record until validated; create a clean tested pipeline separately afterward.

## Verified analysis settings visible in source code

- scRNA-seq QC: `nFeature_RNA > 300`, `nFeature_RNA < 7000`, `percent.mt < 20`
- `NormalizeData(..., normalization.method = "LogNormalize", scale.factor = 10000)`
- `FindVariableFeatures(..., selection.method = "vst", nfeatures = 3000)`
- Ism1 `FeaturePlot` colour limits: `min.cutoff = "q05"`, `max.cutoff = "q99"`; these are **plotting limits, not cell-selection thresholds**

## Before making the GitHub repository public

1. Confirm with the supervisor that these analysis scripts can be disclosed.
2. Review code for hidden tokens, URLs, raw data, private sample identifiers, local paths, and proprietary metadata.
3. Verify individual scripts with actual approved input datasets; capture `sessionInfo()` and version-lock dependencies.
4. Match each figure and statistical cutoff against the final thesis; document which script generated it.
5. Recover or recreate and validate the upstream GSE271812 data-import/object-building workflow; downstream original code is now included.

**Source code only:** No dataset access credentials, sequencing reads, matrices or intermediate `.rds` objects should be committed without explicit approval.
