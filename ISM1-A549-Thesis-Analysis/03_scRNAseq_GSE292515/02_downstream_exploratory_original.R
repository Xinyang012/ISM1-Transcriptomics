# Thesis source-code archive — reviewer draft, NOT validated end-to-end.
# Original analysis logic, thresholds and sample identifiers preserved.
# Replace PATH_TO_* placeholders with permitted local paths before running.
# Check README.md and notes/VALIDATION_CHECKLIST.md before sharing publicly.

##############################################################
## 0. 恢复GSE292515分析环境
##############################################################

rm(list = ls())

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(harmony)
library(DESeq2)

set.seed(1234)

base.dir <- "PATH_TO_GSE292515_DATA"

output.dir <- file.path(
  base.dir,
  "results_new_analysis"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

qc.file <- file.path(
  base.dir,
  "results_Ism1_epithelial",
  "GSE292515_after_QC.rds"
)

file.exists(qc.file)

gse292515 <- readRDS(qc.file)

print(gse292515)

table(
  gse292515$age,
  gse292515$treatment
)



##############################################################
## Normalize
##############################################################

Seurat::DefaultAssay(gse292515) <- "RNA"

if (
  inherits(
    gse292515[["RNA"]],
    "Assay5"
  )
) {
  
  gse292515 <- SeuratObject::JoinLayers(
    gse292515,
    assay = "RNA"
  )
}

gse292515 <- Seurat::NormalizeData(
  gse292515,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

gse292515 <- Seurat::FindVariableFeatures(
  gse292515,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

gse292515 <- Seurat::ScaleData(
  gse292515,
  features = Seurat::VariableFeatures(gse292515),
  verbose = FALSE
)

gse292515 <- Seurat::RunPCA(
  gse292515,
  features = Seurat::VariableFeatures(gse292515),
  npcs = 30,
  verbose = FALSE
)

dims.use <- 1:20




##############################################################
## Harmony
##############################################################

gse292515 <- harmony::RunHarmony(
  object = gse292515,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = dims.use,
  reduction.save = "harmony"
)

gse292515 <- Seurat::FindNeighbors(
  gse292515,
  reduction = "harmony",
  dims = dims.use,
  verbose = FALSE
)

gse292515 <- Seurat::FindClusters(
  gse292515,
  resolution = 0.5,
  random.seed = 1234,
  verbose = FALSE
)

gse292515 <- Seurat::RunUMAP(
  gse292515,
  reduction = "harmony",
  dims = dims.use,
  seed.use = 1234,
  verbose = FALSE
)





sort(
  unique(
    as.numeric(
      as.character(
        gse292515$seurat_clusters
      )
    )
  )
)






##############################################################
## 恢复broad cell annotation
##############################################################

cluster.to.celltype <- c(
  
  "0"  = "Endothelial",
  "1"  = "Fibroblasts",
  "2"  = "Macrophages",
  "3"  = "Epithelial",
  "4"  = "Epithelial",
  "5"  = "B Cells",
  "6"  = "NK Cells",
  "7"  = "Epithelial",
  "8"  = "T Cells",
  "9"  = "Dendritic Cells",
  "10" = "Macrophages",
  "11" = "Ciliated Cells",
  "12" = "Epithelial",
  "13" = "Monocytes",
  "14" = "Fibroblasts",
  "15" = "NK Cells",
  "16" = "Fibroblasts",
  "17" = "Epithelial",
  "18" = "Mixed/Doublet candidate",
  "19" = "Neutrophils",
  "20" = "Fibroblasts",
  "21" = "Epithelial"
)

gse292515$celltype_broad <- unname(
  cluster.to.celltype[
    as.character(
      gse292515$seurat_clusters
    )
  ]
)

table(
  Cluster = gse292515$seurat_clusters,
  CellType = gse292515$celltype_broad
)





saveRDS(
  gse292515,
  file = file.path(
    base.dir,
    "GSE292515_annotated_for_new_analysis.rds"
  )
)





##############################################################
## 1. Ism1-positive fraction
##############################################################

gse.analysis <- subset(
  gse292515,
  subset =
    celltype_broad !=
    "Mixed/Doublet candidate"
)

Seurat::DefaultAssay(gse.analysis) <- "RNA"





ism1.cell.data <- Seurat::FetchData(
  object = gse.analysis,
  vars = c(
    "Ism1",
    "sample",
    "age",
    "treatment",
    "group",
    "seurat_clusters",
    "celltype_broad"
  )
)

ism1.sample.cluster <-
  ism1.cell.data |>
  dplyr::group_by(
    sample,
    age,
    treatment,
    group,
    seurat_clusters,
    celltype_broad
  ) |>
  dplyr::summarise(
    
    n_cells = dplyr::n(),
    
    mean_Ism1 = mean(
      Ism1,
      na.rm = TRUE
    ),
    
    pct_Ism1_positive = mean(
      Ism1 > 0,
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  )


ism1.sample.cluster.use <-
  ism1.sample.cluster |>
  dplyr::filter(
    n_cells >= 20
  )



write.csv(
  ism1.sample.cluster.use,
  file.path(
    output.dir,
    "01_Ism1_positive_fraction_by_mouse_cluster.csv"
  ),
  row.names = FALSE
)




ism1.fraction.group <-
  ism1.sample.cluster.use |>
  dplyr::group_by(
    group,
    seurat_clusters,
    celltype_broad
  ) |>
  dplyr::summarise(
    
    mean_pct_positive =
      mean(
        pct_Ism1_positive,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    
    cluster_label = paste0(
      "C",
      seurat_clusters,
      " | ",
      celltype_broad
    )
  )




p.ism1.fraction <- ggplot2::ggplot(
  ism1.fraction.group,
  ggplot2::aes(
    x = group,
    y = cluster_label,
    fill = mean_pct_positive
  )
) +
  ggplot2::geom_tile() +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "Ism1-positive cells across lung cell clusters",
    x = NULL,
    y = NULL,
    fill = "% Ism1+"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

print(p.ism1.fraction)






ggplot2::ggsave(
  filename = file.path(
    output.dir,
    "01_Ism1_positive_fraction_heatmap.png"
  ),
  plot = p.ism1.fraction,
  width = 9,
  height = 9,
  dpi = 300,
  bg = "white"
)





##############################################################
## 2. Cluster pseudobulk
##############################################################

gse.analysis$cluster_id <- as.character(
  gse.analysis$seurat_clusters
)

gse.analysis$pb_id <- paste0(
  gse.analysis$sample,
  "__C",
  gse.analysis$cluster_id
)



pb.cell.number <-
  gse.analysis@meta.data |>
  dplyr::count(
    pb_id,
    name = "n_cells"
  )


pb.counts <- Seurat::AggregateExpression(
  object = gse.analysis,
  assays = "RNA",
  group.by = "pb_id",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA





pb.meta <-
  gse.analysis@meta.data |>
  dplyr::transmute(
    
    pb_id = pb_id,
    sample = sample,
    age = age,
    treatment = treatment,
    group = group,
    cluster = cluster_id,
    
    celltype_broad =
      as.character(
        celltype_broad
      )
  ) |>
  dplyr::distinct() |>
  dplyr::left_join(
    pb.cell.number,
    by = "pb_id"
  )

print(pb.meta)




clusters.use <- sort(
  unique(
    pb.meta$cluster
  )
)

ism1.pb.list <- list()


##############################################################
## 修正pseudobulk count matrix的列名
##############################################################

colnames(pb.counts) <- gsub(
  "--",
  "__",
  colnames(pb.counts),
  fixed = TRUE
)

## 检查是否修复成功
head(pb.meta$pb_id)
head(colnames(pb.counts))

setdiff(
  pb.meta$pb_id,
  colnames(pb.counts)
)

if (anyDuplicated(colnames(pb.counts)) > 0) {
  stop("修正列名后出现重复的pseudobulk ID。")
}


##############################################################
## Cluster-specific pseudobulk DESeq2
## IAV vs Mock，校正age
##############################################################

ism1.pb.list <- list()

for (current.cluster in clusters.use) {
  
  cat(
    "\nAnalyzing cluster ",
    current.cluster,
    "\n"
  )
  
  current.meta <-
    pb.meta |>
    dplyr::filter(
      cluster == current.cluster,
      n_cells >= 20
    )
  
  
  ############################################################
  ## 基本检查
  ############################################################
  
  if (nrow(current.meta) < 6) {
    cat("Skipped: too few samples.\n")
    next
  }
  
  if (
    dplyr::n_distinct(
      current.meta$treatment
    ) < 2
  ) {
    cat("Skipped: treatment missing.\n")
    next
  }
  
  if (
    dplyr::n_distinct(
      current.meta$age
    ) < 2
  ) {
    cat("Skipped: age group missing.\n")
    next
  }
  
  
  ############################################################
  ## factors
  ############################################################
  
  current.meta$age <- factor(
    current.meta$age,
    levels = c(
      "Young",
      "Old"
    )
  )
  
  current.meta$treatment <- factor(
    current.meta$treatment,
    levels = c(
      "Mock",
      "IAV"
    )
  )
  
  current.meta <- as.data.frame(
    current.meta
  )
  
  rownames(
    current.meta
  ) <- current.meta$pb_id
  
  
  ############################################################
  ## counts
  ############################################################
  
  current.counts <- pb.counts[
    ,
    current.meta$pb_id,
    drop = FALSE
  ]
  
  current.counts <- round(
    as.matrix(
      current.counts
    )
  )
  
  
  ############################################################
  ## design检查
  ############################################################
  
  design.matrix <- stats::model.matrix(
    ~ age + treatment,
    data = current.meta
  )
  
  if (
    qr(design.matrix)$rank <
    ncol(design.matrix)
  ) {
    
    cat(
      "Skipped: design not full rank.\n"
    )
    
    next
  }
  
  
  ############################################################
  ## DESeq2
  ############################################################
  
  dds.temp <-
    DESeq2::DESeqDataSetFromMatrix(
      
      countData = current.counts,
      
      colData = current.meta,
      
      design = ~ age + treatment
    )
  
  
  keep.genes <- rowSums(
    DESeq2::counts(
      dds.temp
    ) >= 5
  ) >= 2
  
  dds.temp <- dds.temp[
    keep.genes,
  ]
  
  
  if (
    !(
      "Ism1" %in%
      rownames(dds.temp)
    )
  ) {
    
    cat(
      "Skipped: Ism1 too low.\n"
    )
    
    next
  }
  
  
  dds.temp <- DESeq2::DESeq(
    dds.temp,
    quiet = TRUE
  )
  
  
  res.temp <- DESeq2::results(
    dds.temp,
    contrast = c(
      "treatment",
      "IAV",
      "Mock"
    )
  )
  
  
  ism1.temp <- as.data.frame(
    res.temp[
      "Ism1",
      ,
      drop = FALSE
    ]
  )
  
  
  ism1.temp$cluster <-
    as.character(
      current.cluster
    )
  
  ism1.temp$celltype_broad <-
    as.character(
      unique(
        current.meta$celltype_broad
      )[1]
    )
  
  ism1.temp$n_samples <-
    nrow(
      current.meta
    )
  
  
  ############################################################
  ## 保存当前cluster结果
  ############################################################
  
  ism1.pb.list[[as.character(current.cluster)]] <- ism1.temp
  
  cat(
    "Finished cluster ",
    current.cluster,
    "\n",
    sep = ""
  )
}




ism1.pb.results <-
  dplyr::bind_rows(
    ism1.pb.list
  ) |>
  dplyr::mutate(
    
    padj_across_clusters =
      stats::p.adjust(
        pvalue,
        method = "BH"
      ),
    
    lower95 =
      log2FoldChange -
      1.96 * lfcSE,
    
    upper95 =
      log2FoldChange +
      1.96 * lfcSE
  ) |>
  dplyr::arrange(
    dplyr::desc(
      log2FoldChange
    )
  )

print(
  ism1.pb.results,
  row.names = FALSE
)



write.csv(
  ism1.pb.results,
  file.path(
    output.dir,
    "02_Ism1_cluster_pseudobulk_IAV_vs_Mock.csv"
  ),
  row.names = FALSE
)





ism1.pb.plot.data <-
  ism1.pb.results |>
  dplyr::mutate(
    
    cluster_label =
      paste0(
        "C",
        cluster,
        " | ",
        celltype_broad
      )
  )




p.ism1.pb <- ggplot2::ggplot(
  ism1.pb.plot.data,
  ggplot2::aes(
    x = log2FoldChange,
    y = reorder(
      cluster_label,
      log2FoldChange
    )
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0
  ) +
  
  ggplot2::geom_point(
    size = 2.8
  ) +
  
  ggplot2::theme_classic() +
  
  ggplot2::labs(
    title = "Cell-type-specific Ism1 response to IAV",
    subtitle = "Sample-level pseudobulk analysis",
    x = "Ism1 log2 fold change: IAV vs Mock",
    y = NULL
  )

print(p.ism1.pb)






ggplot2::ggsave(
  filename = file.path(
    output.dir,
    "02_Ism1_pseudobulk_forest_plot.png"
  ),
  plot = p.ism1.pb,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)


ism1.fraction.treatment <-
  ism1.sample.cluster.use |>
  dplyr::group_by(
    seurat_clusters,
    celltype_broad,
    treatment
  ) |>
  dplyr::summarise(
    
    mean_pct =
      mean(
        pct_Ism1_positive,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = mean_pct,
    names_prefix = "pct_"
  ) |>
  dplyr::mutate(
    
    delta_pct =
      pct_IAV -
      pct_Mock
  )








cluster.priority <-
  ism1.pb.results |>
  dplyr::left_join(
    ism1.fraction.treatment,
    by = c(
      "cluster" =
        "seurat_clusters",
      
      "celltype_broad" =
        "celltype_broad"
    )
  ) |>
  dplyr::select(
    cluster,
    celltype_broad,
    log2FoldChange,
    pvalue,
    padj,
    padj_across_clusters,
    pct_Mock,
    pct_IAV,
    delta_pct
  ) |>
  dplyr::arrange(
    padj_across_clusters,
    dplyr::desc(
      abs(
        log2FoldChange
      )
    )
  )

print(
  cluster.priority,
  row.names = FALSE
)






write.csv(
  cluster.priority,
  file.path(
    output.dir,
    "03_Ism1_key_cluster_priority.csv"
  ),
  row.names = FALSE
)














##############################################################
## 4A. 定义重点Ism1-related clusters
##############################################################

key.clusters <- c(
  "1",   # Fibroblasts
  "9",   # Dendritic Cells
  "7"    # Epithelial - exploratory
)

cluster.label.map <- c(
  "1" = "C1 | Fibroblasts",
  "9" = "C9 | Dendritic Cells",
  "7" = "C7 | Epithelial"
)

cat(
  "Key clusters:\n"
)

print(
  cluster.label.map
)
  
  
  
  
  
  
  
##############################################################
## 4B. Core antiviral ISG signature
##############################################################

core.antiviral.genes <- c(
  "Isg15",
  "Ifit1",
  "Ifit2",
  "Ifit3",
  "Mx1",
  "Mx2",
  "Oas1a",
  "Oas2",
  "Oas3",
  "Rsad2",
  "Usp18",
  "Irf7",
  "Ddx58"
)

core.antiviral.present <- intersect(
  core.antiviral.genes,
  rownames(gse.analysis)
)

cat(
  "Genes available for antiviral score:\n"
)

print(
  core.antiviral.present
)

cat(
  "\nNumber of genes:",
  length(core.antiviral.present),
  "\n"
)

if (length(core.antiviral.present) < 5) {
  
  stop(
    "Too few antiviral genes are present."
  )
}  
  
  
  


##############################################################
## 删除旧Antiviral_ISG score
##############################################################

if (
  "Antiviral_ISG1" %in%
  colnames(gse.analysis@meta.data)
) {
  
  gse.analysis$Antiviral_ISG1 <- NULL
}


##############################################################
## 计算antiviral ISG module score
##############################################################

Seurat::DefaultAssay(
  gse.analysis
) <- "RNA"

gse.analysis <- Seurat::AddModuleScore(
  object = gse.analysis,
  features = list(
    core.antiviral.present
  ),
  assay = "RNA",
  name = "Antiviral_ISG",
  seed = 1234
)

summary(
  gse.analysis$Antiviral_ISG1
)




saveRDS(
  gse.analysis,
  file = file.path(
    output.dir,
    "GSE292515_with_antiviral_ISG_score.rds"
  )
)



##############################################################
## 4C. 提取重点clusters
##############################################################

key.cells <- rownames(
  gse.analysis@meta.data
)[
  as.character(
    gse.analysis$seurat_clusters
  ) %in%
    key.clusters
]

gse.key <- subset(
  gse.analysis,
  cells = key.cells
)

table(
  Cluster = gse.key$seurat_clusters,
  Treatment = gse.key$treatment
)



gse.key$cluster_label <- unname(
  cluster.label.map[
    as.character(
      gse.key$seurat_clusters
    )
  ]
)

gse.key$cluster_label <- factor(
  gse.key$cluster_label,
  levels = c(
    "C1 | Fibroblasts",
    "C9 | Dendritic Cells",
    "C7 | Epithelial"
  )
)

table(
  gse.key$cluster_label
)





##############################################################
## 4D. Mouse-level antiviral score
##############################################################

antiviral.cell.data <- Seurat::FetchData(
  object = gse.key,
  vars = c(
    "Antiviral_ISG1",
    "sample",
    "age",
    "treatment",
    "group",
    "seurat_clusters",
    "cluster_label"
  )
)




antiviral.sample.cluster <-
  antiviral.cell.data |>
  dplyr::group_by(
    sample,
    age,
    treatment,
    group,
    seurat_clusters,
    cluster_label
  ) |>
  dplyr::summarise(
    
    n_cells = dplyr::n(),
    
    mean_antiviral = mean(
      Antiviral_ISG1,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_cells >= 20
  )

print(
  antiviral.sample.cluster
)




write.csv(
  antiviral.sample.cluster,
  file = file.path(
    output.dir,
    "04_Key_clusters_antiviral_score_by_mouse.csv"
  ),
  row.names = FALSE
)





p.antiviral.treatment <- ggplot2::ggplot(
  antiviral.sample.cluster,
  ggplot2::aes(
    x = treatment,
    y = mean_antiviral,
    shape = age
  )
) +
  
  ggplot2::geom_point(
    size = 3,
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  ggplot2::stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45
  ) +
  
  ggplot2::facet_wrap(
    ~ cluster_label,
    scales = "free_y"
  ) +
  
  ggplot2::theme_classic() +
  
  ggplot2::labs(
    title =
      "Antiviral response in Ism1-related cell populations",
    
    x = NULL,
    
    y =
      "Antiviral ISG module score",
    
    shape =
      "Age"
  )

print(
  p.antiviral.treatment
)



ggplot2::ggsave(
  filename = file.path(
    output.dir,
    "04_Antiviral_score_Mock_vs_IAV_key_clusters.png"
  ),
  plot = p.antiviral.treatment,
  width = 11,
  height = 5.5,
  dpi = 300,
  bg = "white"
)










##############################################################
## 4E. Ism1-positive vs Ism1-negative
##############################################################

ism1.ifn.cells <- Seurat::FetchData(
  object = gse.key,
  vars = c(
    "Ism1",
    "Antiviral_ISG1",
    "sample",
    "age",
    "treatment",
    "group",
    "seurat_clusters",
    "cluster_label"
  )
)






ism1.ifn.cells <-
  ism1.ifn.cells |>
  dplyr::mutate(
    
    Ism1_status =
      dplyr::if_else(
        Ism1 > 0,
        "Ism1_positive",
        "Ism1_negative"
      )
  )




table(
  Cluster = ism1.ifn.cells$seurat_clusters,
  Ism1 = ism1.ifn.cells$Ism1_status
)



##############################################################
## 4F. Mouse-level Ism1+ vs Ism1-
##############################################################

min.cells.per.status <- 5

ism1.ifn.sample <-
  ism1.ifn.cells |>
  dplyr::group_by(
    sample,
    age,
    treatment,
    group,
    seurat_clusters,
    cluster_label,
    Ism1_status
  ) |>
  dplyr::summarise(
    
    n_cells =
      dplyr::n(),
    
    mean_antiviral =
      mean(
        Antiviral_ISG1,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_cells >=
      min.cells.per.status
  )

print(
  ism1.ifn.sample
)



ism1.ifn.delta <-
  ism1.ifn.sample |>
  tidyr::pivot_wider(
    
    names_from =
      Ism1_status,
    
    values_from = c(
      n_cells,
      mean_antiviral
    )
  ) |>
  dplyr::filter(
    
    !is.na(
      mean_antiviral_Ism1_positive
    ),
    
    !is.na(
      mean_antiviral_Ism1_negative
    )
  ) |>
  dplyr::mutate(
    
    delta_antiviral =
      mean_antiviral_Ism1_positive -
      mean_antiviral_Ism1_negative
  )

print(
  ism1.ifn.delta
)


ism1.ifn.delta <-
  ism1.ifn.sample |>
  tidyr::pivot_wider(
    
    names_from =
      Ism1_status,
    
    values_from = c(
      n_cells,
      mean_antiviral
    )
  ) |>
  dplyr::filter(
    
    !is.na(
      mean_antiviral_Ism1_positive
    ),
    
    !is.na(
      mean_antiviral_Ism1_negative
    )
  ) |>
  dplyr::mutate(
    
    delta_antiviral =
      mean_antiviral_Ism1_positive -
      mean_antiviral_Ism1_negative
  )

print(
  ism1.ifn.delta
)




write.csv(
  ism1.ifn.delta,
  file = file.path(
    output.dir,
    "05_Ism1_positive_vs_negative_antiviral_delta.csv"
  ),
  row.names = FALSE
)




##############################################################
## 4I. Type III IFN axis
##############################################################

type3.ifn.genes <- c(
  "Ifnl2",
  "Ifnl3",
  "Ifnlr1",
  "Il10rb"
)

type3.ifn.present <- intersect(
  type3.ifn.genes,
  rownames(gse.key)
)

cat(
  "Type III IFN genes available:\n"
)

print(
  type3.ifn.present
)


gse.key$cluster_treatment <- paste(
  gse.key$cluster_label,
  gse.key$treatment,
  sep = " | "
)







if (length(type3.ifn.present) > 0) {
  
  p.type3.ifn <- Seurat::DotPlot(
    object = gse.key,
    features = type3.ifn.present,
    group.by = "cluster_treatment",
    dot.scale = 7
  ) +
    
    Seurat::RotatedAxis() +
    
    ggplot2::theme_classic() +
    
    ggplot2::labs(
      
      title =
        "Type III IFN axis in Ism1-related cell populations",
      
      x = NULL,
      
      y = NULL
    )
  
  
  print(
    p.type3.ifn
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      output.dir,
      "06_TypeIII_IFN_axis_key_clusters.png"
    ),
    plot = p.type3.ifn,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}




##############################################################
## 5A. 修正pseudobulk IDs
##############################################################

colnames(pb.counts) <- gsub(
  "--",
  "__",
  colnames(pb.counts),
  fixed = TRUE
)

missing.pb <- setdiff(
  pb.meta$pb_id,
  colnames(pb.counts)
)

print(
  missing.pb
)






##############################################################
## 5B. Age × treatment interaction
##############################################################

run.age.interaction <- function(
    cluster.id,
    cluster.label
) {
  
  
  cat(
    "\n====================================\n",
    "Analyzing ",
    cluster.label,
    "\n",
    sep = ""
  )
  
  
  ############################################################
  ## metadata
  ############################################################
  
  current.meta <-
    pb.meta |>
    dplyr::filter(
      as.character(cluster) == cluster.id,
      n_cells >= 20
    )
  
  
  current.meta <- current.meta |>
    dplyr::mutate(
      
      age = factor(
        age,
        levels = c(
          "Young",
          "Old"
        )
      ),
      
      treatment = factor(
        treatment,
        levels = c(
          "Mock",
          "IAV"
        )
      )
    )
  
  
  print(
    table(
      current.meta$age,
      current.meta$treatment
    )
  )
  
  
  ############################################################
  ## 四组都必须存在
  ############################################################
  
  group.table <- table(
    current.meta$age,
    current.meta$treatment
  )
  
  
  if (any(group.table == 0)) {
    
    stop(
      paste0(
        cluster.label,
        ": one age-treatment group is missing."
      )
    )
  }
  
  
  current.meta <- as.data.frame(
    current.meta
  )
  
  
  rownames(
    current.meta
  ) <- current.meta$pb_id
  
  
  ############################################################
  ## counts
  ############################################################
  
  current.counts <- pb.counts[
    ,
    current.meta$pb_id,
    drop = FALSE
  ]
  
  
  current.counts <- round(
    as.matrix(
      current.counts
    )
  )
  
  
  ############################################################
  ## DESeq2
  ############################################################
  
  dds.age <- DESeq2::DESeqDataSetFromMatrix(
    
    countData =
      current.counts,
    
    colData =
      current.meta,
    
    design =
      ~ age * treatment
  )
  
  
  keep.genes <- rowSums(
    DESeq2::counts(
      dds.age
    ) >= 5
  ) >= 2
  
  
  dds.age <- dds.age[
    keep.genes,
  ]
  
  
  if (
    !(
      "Ism1" %in%
      rownames(dds.age)
    )
  ) {
    
    stop(
      paste0(
        cluster.label,
        ": Ism1 was filtered out."
      )
    )
  }
  
  
  dds.age <- DESeq2::DESeq(
    dds.age,
    quiet = TRUE
  )
  
  
  ############################################################
  ## coefficient names
  ############################################################
  
  coefficient.names <-
    DESeq2::resultsNames(
      dds.age
    )
  
  
  print(
    coefficient.names
  )
  
  
  interaction.name <-
    coefficient.names[
      grepl(
        "ageOld.*treatmentIAV|treatmentIAV.*ageOld",
        coefficient.names
      )
    ]
  
  
  if (length(interaction.name) != 1) {
    
    stop(
      paste0(
        cluster.label,
        ": interaction coefficient cannot be uniquely identified."
      )
    )
  }
  
  
  ############################################################
  ## interaction result
  ############################################################
  
  res.interaction <- DESeq2::results(
    dds.age,
    name = interaction.name
  )
  
  
  ism1.interaction <- as.data.frame(
    res.interaction[
      "Ism1",
      ,
      drop = FALSE
    ]
  )
  
  
  ism1.interaction$cluster <-
    cluster.id
  
  ism1.interaction$cluster_label <-
    cluster.label
  
  
  ############################################################
  ## normalized counts for plotting
  ############################################################
  
  normalized.counts <-
    DESeq2::counts(
      dds.age,
      normalized = TRUE
    )
  
  
  plot.data <- data.frame(
    
    sample =
      colnames(
        normalized.counts
      ),
    
    Ism1_normalized =
      as.numeric(
        normalized.counts[
          "Ism1",
        ]
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  plot.data$Ism1_log2 <-
    log2(
      plot.data$Ism1_normalized +
        1
    )
  
  
  meta.for.plot <- data.frame(
    
    sample =
      rownames(
        current.meta
      ),
    
    age =
      current.meta$age,
    
    treatment =
      current.meta$treatment,
    
    stringsAsFactors =
      FALSE
  )
  
  
  plot.data <-
    plot.data |>
    dplyr::left_join(
      meta.for.plot,
      by = "sample"
    )
  
  
  plot.data$cluster <-
    cluster.id
  
  plot.data$cluster_label <-
    cluster.label
  
  
  return(
    list(
      
      summary =
        ism1.interaction,
      
      plot_data =
        plot.data,
      
      dds =
        dds.age
    )
  )
}





##############################################################
## 5C. Run age interactions
##############################################################

age.result.C1 <- run.age.interaction(
  cluster.id = "1",
  cluster.label = "C1 | Fibroblasts"
)

age.result.C9 <- run.age.interaction(
  cluster.id = "9",
  cluster.label = "C9 | Dendritic Cells"
)

age.result.C7 <- run.age.interaction(
  cluster.id = "7",
  cluster.label = "C7 | Epithelial"
)





current.meta <-
  pb.meta |>
  dplyr::filter(
    as.character(cluster) == "7",
    n_cells >= 20
  )

current.counts <- pb.counts[
  ,
  current.meta$pb_id,
  drop = FALSE
]

current.counts <- round(
  as.matrix(current.counts)
)

data.frame(
  sample = current.meta$sample,
  age = current.meta$age,
  treatment = current.meta$treatment,
  Ism1_count = as.numeric(
    current.counts[
      "Ism1",
    ]
  )
)







##############################################################
## 5D. Broad epithelial pseudobulk
## 合并所有Epithelial clusters
##############################################################

epithelial.clusters <- c(
  "3",
  "4",
  "7",
  "12",
  "17",
  "21"
)

epithelial.cells <- rownames(
  gse.analysis@meta.data
)[
  as.character(
    gse.analysis$seurat_clusters
  ) %in%
    epithelial.clusters
]

gse.epithelial <- subset(
  gse.analysis,
  cells = epithelial.cells
)

cat(
  "Broad epithelial cells:",
  ncol(gse.epithelial),
  "\n"
)

table(
  gse.epithelial$sample
)



##############################################################
## 每只mouse聚合Epithelial raw counts
##############################################################

epithelial.aggregate <- Seurat::AggregateExpression(
  object = gse.epithelial,
  assays = "RNA",
  group.by = "sample",
  return.seurat = FALSE,
  verbose = FALSE
)

epithelial.pb.counts <- epithelial.aggregate$RNA

dim(
  epithelial.pb.counts
)

colnames(
  epithelial.pb.counts
)




##############################################################
## Broad epithelial metadata
##############################################################

epithelial.meta <-
  gse.epithelial@meta.data |>
  dplyr::transmute(
    sample = as.character(sample),
    age = as.character(age),
    treatment = as.character(treatment),
    group = as.character(group)
  ) |>
  dplyr::distinct()



epithelial.cell.number <-
  gse.epithelial@meta.data |>
  dplyr::count(
    sample,
    name = "n_cells"
  )

epithelial.meta <-
  epithelial.meta |>
  dplyr::left_join(
    epithelial.cell.number,
    by = "sample"
  )






epithelial.meta <-
  epithelial.meta[
    match(
      colnames(epithelial.pb.counts),
      epithelial.meta$sample
    ),
    ,
    drop = FALSE
  ]

rownames(
  epithelial.meta
) <- epithelial.meta$sample

stopifnot(
  identical(
    colnames(epithelial.pb.counts),
    rownames(epithelial.meta)
  )
)

print(
  epithelial.meta
)







##############################################################
## 检查Broad epithelial Ism1 raw counts
##############################################################

epithelial.ism1.raw <- data.frame(
  
  sample =
    epithelial.meta$sample,
  
  age =
    epithelial.meta$age,
  
  treatment =
    epithelial.meta$treatment,
  
  n_cells =
    epithelial.meta$n_cells,
  
  Ism1_count =
    as.numeric(
      epithelial.pb.counts[
        "Ism1",
      ]
    ),
  
  stringsAsFactors =
    FALSE
)

print(
  epithelial.ism1.raw
)


##############################################################
## Broad epithelial age × treatment interaction
##############################################################

epithelial.meta$age <- factor(
  epithelial.meta$age,
  levels = c(
    "Young",
    "Old"
  )
)

epithelial.meta$treatment <- factor(
  epithelial.meta$treatment,
  levels = c(
    "Mock",
    "IAV"
  )
)

print(
  table(
    epithelial.meta$age,
    epithelial.meta$treatment
  )
)




dds.epithelial <- DESeq2::DESeqDataSetFromMatrix(
  
  countData = round(
    as.matrix(
      epithelial.pb.counts
    )
  ),
  
  colData = epithelial.meta,
  
  design = ~ age * treatment
)




keep.epithelial <- rowSums(
  DESeq2::counts(
    dds.epithelial
  ) >= 5
) >= 2

if (
  "Ism1" %in%
  rownames(dds.epithelial)
) {
  
  keep.epithelial["Ism1"] <- TRUE
}

dds.epithelial <-
  dds.epithelial[
    keep.epithelial,
  ]

"Ism1" %in%
  rownames(dds.epithelial)



dds.epithelial <- DESeq2::DESeq(
  dds.epithelial,
  quiet = TRUE
)

DESeq2::resultsNames(
  dds.epithelial
)


epithelial.coef.names <-
  DESeq2::resultsNames(
    dds.epithelial
  )

epithelial.interaction.name <-
  epithelial.coef.names[
    grepl(
      "ageOld.*treatmentIAV|treatmentIAV.*ageOld",
      epithelial.coef.names
    )
  ]

print(
  epithelial.interaction.name
)

if (
  length(
    epithelial.interaction.name
  ) != 1
) {
  
  stop(
    "无法唯一识别Broad Epithelial的age × treatment interaction coefficient。"
  )
}



res.epithelial.interaction <-
  DESeq2::results(
    dds.epithelial,
    name =
      epithelial.interaction.name
  )

ism1.epithelial.interaction <-
  as.data.frame(
    res.epithelial.interaction[
      "Ism1",
      ,
      drop = FALSE
    ]
  )

ism1.epithelial.interaction$cluster <-
  "Broad_Epithelial"

ism1.epithelial.interaction$cluster_label <-
  "Broad Epithelial"

print(
  ism1.epithelial.interaction
)





##############################################################
## 4. Spearman correlation
##############################################################

cor.test.ifn.271812 <-
  cor.test(
    gsva.plot.data.271812$Ism1_VST,
    gsva.plot.data.271812$GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION,
    method = "spearman",
    exact = FALSE
  )

print(
  cor.test.ifn.271812
)













##############################################################
## Young-only / Old-only cluster-specific pseudobulk
## Ism1: IAV vs Mock
##############################################################
##############################################################
## Young-only / Old-only cluster pseudobulk Ism1 response
##############################################################

library(DESeq2)
library(dplyr)
library(ggplot2)

run.age.specific.ism1 <- function(age.use) {
  
  result.list <- list()
  
  for (current.cluster in clusters.use) {
    
    cat(
      "\n",
      age.use,
      " | cluster ",
      current.cluster,
      "\n",
      sep = ""
    )
    
    ############################################################
    ## Metadata
    ############################################################
    
    current.meta <- pb.meta |>
      dplyr::filter(
        cluster == current.cluster,
        age == age.use,
        n_cells >= 20
      )
    
    
    # 至少需要两个 treatment group
    if (
      dplyr::n_distinct(current.meta$treatment) < 2
    ) {
      cat("Skipped: treatment group missing.\n")
      next
    }
    
    
    # 每组至少建议有2个 biological samples
    sample.table <- table(current.meta$treatment)
    
    if (
      any(sample.table[c("Mock", "IAV")] < 2)
    ) {
      cat("Skipped: too few biological replicates.\n")
      next
    }
    
    
    current.meta$treatment <- factor(
      current.meta$treatment,
      levels = c("Mock", "IAV")
    )
    
    current.meta <- as.data.frame(current.meta)
    
    rownames(current.meta) <- current.meta$pb_id
    
    
    ############################################################
    ## Counts
    ############################################################
    
    current.counts <- pb.counts[
      ,
      current.meta$pb_id,
      drop = FALSE
    ]
    
    current.counts <- round(
      as.matrix(current.counts)
    )
    
    
    ############################################################
    ## DESeq2
    ############################################################
    
    dds.temp <- DESeq2::DESeqDataSetFromMatrix(
      countData = current.counts,
      colData = current.meta,
      design = ~ treatment
    )
    
    
    keep.genes <- rowSums(
      DESeq2::counts(dds.temp) >= 5
    ) >= 2
    
    dds.temp <- dds.temp[
      keep.genes,
    ]
    
    
    if (!("Ism1" %in% rownames(dds.temp))) {
      cat("Skipped: Ism1 too low.\n")
      next
    }
    
    
    dds.temp <- DESeq2::DESeq(
      dds.temp,
      quiet = TRUE
    )
    
    
    res.temp <- DESeq2::results(
      dds.temp,
      contrast = c(
        "treatment",
        "IAV",
        "Mock"
      )
    )
    
    
    ism1.temp <- as.data.frame(
      res.temp[
        "Ism1",
        ,
        drop = FALSE
      ]
    )
    
    
    ############################################################
    ## Information
    ############################################################
    
    ism1.temp$cluster <- as.character(
      current.cluster
    )
    
    ism1.temp$celltype_broad <- as.character(
      unique(current.meta$celltype_broad)[1]
    )
    
    ism1.temp$age <- age.use
    
    ism1.temp$n_Mock <- sum(
      current.meta$treatment == "Mock"
    )
    
    ism1.temp$n_IAV <- sum(
      current.meta$treatment == "IAV"
    )
    
    
    result.list[[as.character(current.cluster)]] <-
      ism1.temp
  }
  
  
  ############################################################
  ## Combine
  ############################################################
  
  result.df <- dplyr::bind_rows(
    result.list
  )
  
  
  result.df <- result.df |>
    dplyr::mutate(
      
      # across-cluster BH adjustment
      padj_across_clusters =
        stats::p.adjust(
          pvalue,
          method = "BH"
        ),
      
      lower95 =
        log2FoldChange -
        1.96 * lfcSE,
      
      upper95 =
        log2FoldChange +
        1.96 * lfcSE,
      
      cluster_label =
        paste0(
          "C",
          cluster,
          " | ",
          celltype_broad
        )
    ) |>
    dplyr::arrange(
      dplyr::desc(log2FoldChange)
    )
  
  
  return(result.df)
}



##############################################################
## Young-only
##############################################################

ism1.pb.young <- run.age.specific.ism1(
  age.use = "Young"
)


##############################################################
## Old-only
##############################################################

ism1.pb.old <- run.age.specific.ism1(
  age.use = "Old"
)







ism1.pb.young <- ism1.pb.young |>
  dplyr::mutate(
    significance = dplyr::case_when(
      padj_across_clusters < 0.05 &
        log2FoldChange > 0 ~ "Upregulated",
      
      padj_across_clusters < 0.05 &
        log2FoldChange < 0 ~ "Downregulated",
      
      TRUE ~ "Not significant"
    )
  )

ism1.pb.old <- ism1.pb.old |>
  dplyr::mutate(
    significance = dplyr::case_when(
      padj_across_clusters < 0.05 &
        log2FoldChange > 0 ~ "Upregulated",
      
      padj_across_clusters < 0.05 &
        log2FoldChange < 0 ~ "Downregulated",
      
      TRUE ~ "Not significant"
    )
  )




p.ism1.young <- ggplot(
  ism1.pb.young,
  aes(
    x = log2FoldChange,
    y = reorder(
      cluster_label,
      log2FoldChange
    )
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0,
    linewidth = 0.7,
    color = "grey35"
  ) +
  
  geom_point(
    aes(
      fill = significance
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.45,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Upregulated" = "#D95F59",
      "Downregulated" = "#4C78A8",
      "Not significant" = "white"
    )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 8)
    ),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "grey35"
    ),
    plot.margin = margin(
      10, 15, 10, 10
    )
  ) +
  
  labs(
    title = "Young mice",
    subtitle =
      "Cell-type-specific Ism1 response to IAV",
    x = expression(
      Ism1~log[2]~fold~change~"(IAV vs Mock)"
    ),
    y = NULL
  )

p.ism1.young



p.ism1.old <- ggplot(
  ism1.pb.old,
  aes(
    x = log2FoldChange,
    y = reorder(
      cluster_label,
      log2FoldChange
    )
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0,
    linewidth = 0.7,
    color = "grey35"
  ) +
  
  geom_point(
    aes(
      fill = significance
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.45,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Upregulated" = "#D95F59",
      "Downregulated" = "#4C78A8",
      "Not significant" = "white"
    )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 8)
    ),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "grey35"
    ),
    plot.margin = margin(
      10, 15, 10, 10
    )
  ) +
  
  labs(
    title = "Old mice",
    subtitle =
      "Cell-type-specific Ism1 response to IAV",
    x = expression(
      Ism1~log[2]~fold~change~"(IAV vs Mock)"
    ),
    y = NULL
  )

p.ism1.old







library(dplyr)
library(ggplot2)

ism1.pb.pretty <- ism1.pb.results |>
  dplyr::mutate(
    
    cluster_label = paste0(
      "C",
      cluster,
      " | ",
      celltype_broad
    ),
    
    significance = dplyr::case_when(
      
      padj_across_clusters < 0.05 &
        log2FoldChange > 0 ~
        "Significantly upregulated",
      
      padj_across_clusters < 0.05 &
        log2FoldChange < 0 ~
        "Significantly downregulated",
      
      TRUE ~
        "Not significant"
    ),
    
    sig_label = dplyr::case_when(
      padj_across_clusters < 0.001 ~ "***",
      padj_across_clusters < 0.01  ~ "**",
      padj_across_clusters < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) |>
  
  dplyr::arrange(
    log2FoldChange
  )

ism1.pb.pretty$cluster_label <- factor(
  ism1.pb.pretty$cluster_label,
  levels = ism1.pb.pretty$cluster_label
)




p.ism1.pb.pretty <- ggplot(
  ism1.pb.pretty,
  aes(
    x = log2FoldChange,
    y = cluster_label
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey45"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0,
    linewidth = 0.75,
    color = "grey35"
  ) +
  
  geom_point(
    aes(
      fill = significance
    ),
    shape = 21,
    size = 4.2,
    stroke = 0.55,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Significantly upregulated" = "#D95F59",
      "Significantly downregulated" = "#4C78A8",
      "Not significant" = "#F2F2F2"
    )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    axis.text.y = element_text(
      size = 11.5,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 10)
    ),
    
    axis.line = element_line(
      linewidth = 0.7
    ),
    
    legend.position = "top",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 10.5
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 11.5,
      color = "grey35",
      margin = margin(b = 12)
    ),
    
    plot.margin = margin(
      15, 25, 15, 15
    )
  ) +
  
  labs(
    title =
      "Cell-type-specific Ism1 response to IAV",
    
    subtitle =
      "Sample-level pseudobulk DESeq2 analysis adjusted for age",
    
    x =
      expression(
        Ism1~log[2]~fold~change~"(IAV vs Mock)"
      ),
    
    y = NULL
  )

p.ism1.pb.pretty




ism1.pb.pretty <- ism1.pb.pretty |>
  dplyr::mutate(
    label_x = ifelse(
      log2FoldChange >= 0,
      upper95 + 0.12,
      lower95 - 0.12
    )
  )




p.ism1.pb.pretty <-
  p.ism1.pb.pretty +
  
  geom_text(
    data = subset(
      ism1.pb.pretty,
      sig_label != ""
    ),
    aes(
      x = label_x,
      y = cluster_label,
      label = sig_label
    ),
    size = 4.5,
    fontface = "bold",
    inherit.aes = FALSE
  )

p.ism1.pb.pretty









ggsave(
  filename = file.path(
    output.dir,
    "02_Ism1_pseudobulk_forest_plot_pretty.png"
  ),
  plot = p.ism1.pb.pretty,
  width = 8.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)
