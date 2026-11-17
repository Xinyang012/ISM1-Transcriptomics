# Thesis source-code archive — reviewer draft, NOT validated end-to-end.
# Original analysis logic, thresholds and sample identifiers preserved.
# Replace PATH_TO_* placeholders with permitted local paths before running.
# Check README.md and notes/VALIDATION_CHECKLIST.md before sharing publicly.

## GSE292515：下载并检查文件
##############################################################

rm(list = ls())

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

setwd(
  "PATH_TO_SCRNASEQ_WORKDIR"
)

dir.create(
  "GSE292515_MOUSE",
  showWarnings = FALSE
)

setwd("GSE292515_MOUSE")

##############################################################
## 下载GEO补充文件
##############################################################

getGEOSuppFiles(
  GEO = "GSE292515",
  makeDirectory = TRUE,
  baseDir = "."
)

##############################################################
## 查看下载到的全部文件
##############################################################

downloaded.files <- list.files(
  "GSE292515",
  recursive = TRUE,
  full.names = TRUE
)

file.table <- data.frame(
  file = basename(downloaded.files),
  path = downloaded.files,
  size_MB = round(
    file.info(downloaded.files)$size /
      1024^2,
    2
  ),
  stringsAsFactors = FALSE
)

print(
  file.table,
  row.names = FALSE
)

write.csv(
  file.table,
  "GSE292515_downloaded_files.csv",
  row.names = FALSE
)



##############################################################
## GSE292515：检查H5文件
##############################################################

rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(1234)

setwd(
  "PATH_TO_SCRNASEQ_WORKDIR/GSE292515"
)

raw.dir <- "RAW"

h5.files <- list.files(
  path = raw.dir,
  pattern = "_filtered_feature_bc_matrix\\.h5$",
  full.names = TRUE
)

h5.files <- sort(h5.files)

cat("找到的H5文件数量：", length(h5.files), "\n")

print(
  basename(h5.files)
)

install.packages("hdf5r")
library(hdf5r)

packageVersion("hdf5r")

##############################################################
## 测试读取第一个H5文件
##############################################################

test.counts <- Read10X_h5(
  filename = h5.files[1],
  use.names = TRUE,
  unique.features = TRUE
)

class(test.counts)


dim(test.counts)

head(rownames(test.counts))
head(colnames(test.counts))

test.obj <- CreateSeuratObject(
  counts = test.counts,
  project = "test_sample",
  min.cells = 3,
  min.features = 200
)

print(test.obj)



h5.files <- list.files(
  raw.dir,
  pattern = "_filtered_feature_bc_matrix\\.h5$",
  full.names = TRUE
)

h5.files <- sort(h5.files)

length(h5.files)
basename(h5.files)


##############################################################
## 提取样本名称
##############################################################

sample.names <- sub(
  "_filtered_feature_bc_matrix\\.h5$",
  "",
  basename(h5.files)
)

sample.names




##############################################################
## 建立GSM和实验分组对应表
##############################################################

##############################################################
## 建立GSM和实验分组对应表
##############################################################

sample.meta <- data.frame(
  
  GSM = c(
    "GSM8860478",
    "GSM8860479",
    "GSM8860480",
    "GSM8860481",
    "GSM8860482",
    "GSM8860483",
    "GSM8860484",
    "GSM8860485",
    "GSM8860486",
    "GSM8860487"
  ),
  
  age = c(
    "Young",
    "Young",
    "Young",
    "Young",
    "Young",
    "Old",
    "Old",
    "Old",
    "Old",
    "Old"
  ),
  
  treatment = c(
    "Mock",
    "Mock",
    "IAV",
    "IAV",
    "IAV",
    "Mock",
    "Mock",
    "IAV",
    "IAV",
    "IAV"
  ),
  
  replicate = c(
    1, 2,
    1, 2, 3,
    1, 2,
    1, 2, 3
  ),
  
  stringsAsFactors = FALSE
)

sample.meta$group <- paste(
  sample.meta$age,
  sample.meta$treatment,
  sep = "_"
)

print(sample.meta)




##############################################################
## 从文件名提取GSM编号
##############################################################

gsm.names <- sub(
  "^(GSM[0-9]+)_.*$",
  "\\1",
  sample.names
)

gsm.names

## 检查每个文件是否都能在sample.meta中找到
setdiff(
  gsm.names,
  sample.meta$GSM
)



##############################################################
## 读取全部样本
##############################################################

sample.list <- list()

for (i in seq_along(h5.files)) {
  
  sample <- sample.names[i]
  gsm <- gsm.names[i]
  
  cat(
    "[",
    i,
    "/",
    length(h5.files),
    "] Reading: ",
    sample,
    "\n",
    sep = ""
  )
  
  counts <- Read10X_h5(
    filename = h5.files[i],
    use.names = TRUE,
    unique.features = TRUE
  )
  
  ## 正常情况下你的文件返回dgCMatrix
  ## 若某个文件意外返回list，只保留Gene Expression
  if (is.list(counts)) {
    
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = gsm,
    min.cells = 3,
    min.features = 200
  )
  
  ## 给细胞barcode加样本前缀
  obj <- RenameCells(
    obj,
    add.cell.id = gsm
  )
  
  ## 找到当前样本对应的元数据行
  meta.row <- sample.meta[
    sample.meta$GSM == gsm,
    ,
    drop = FALSE
  ]
  
  if (nrow(meta.row) != 1) {
    stop(
      "无法唯一匹配样本元数据：",
      gsm
    )
  }
  
  ## 添加元数据
  obj$sample <- gsm
  obj$sample_file <- sample
  obj$age <- meta.row$age
  obj$treatment <- meta.row$treatment
  obj$replicate <- meta.row$replicate
  obj$group <- meta.row$group
  
  sample.list[[gsm]] <- obj
  
  cat(
    "Genes:",
    nrow(obj),
    " Cells:",
    ncol(obj),
    " Group:",
    meta.row$group,
    "\n"
  )
}

##############################################################
## 检查各样本细胞数
##############################################################

sample.cell.number <- data.frame(
  sample = names(sample.list),
  n_cells = sapply(sample.list, ncol),
  stringsAsFactors = FALSE
)

sample.cell.number <- left_join(
  sample.cell.number,
  sample.meta,
  by = c("sample" = "GSM")
)

print(sample.cell.number)

sample.list[[1]]

head(
  sample.list[[1]]@meta.data
)




##############################################################
## 合并全部样本
##############################################################

gse292515 <- merge(
  x = sample.list[[1]],
  y = sample.list[-1],
  project = "GSE292515"
)

print(gse292515)

dim(gse292515)



##############################################################
## 检查样本分组是否正确
##############################################################

table(
  gse292515$sample,
  gse292515$group
)

table(
  Age = gse292515$age,
  Treatment = gse292515$treatment
)


anyDuplicated(
  colnames(gse292515)
)


saveRDS(
  gse292515,
  file = "GSE292515_raw_merged.rds"
)



##############################################################
## 计算线粒体和核糖体基因比例
##############################################################

gse292515[["percent.mt"]] <- PercentageFeatureSet(
  gse292515,
  pattern = "^mt-"
)

gse292515[["percent.ribo"]] <- PercentageFeatureSet(
  gse292515,
  pattern = "^Rp[sl]"
)





summary(gse292515$percent.mt)
summary(gse292515$percent.ribo)

head(
  gse292515@meta.data[
    ,
    c(
      "sample",
      "age",
      "treatment",
      "group",
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.ribo"
    )
  ]
)


head(
  rownames(gse292515),
  30
)


grep(
  "^Rp[sl]",
  rownames(gse292515),
  value = TRUE
) |>
  head(20)




p.qc.before <- VlnPlot(
  gse292515,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.ribo"
  ),
  group.by = "sample",
  pt.size = 0,
  ncol = 4,
  layer = "counts"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.title.x = element_blank()
  )

print(p.qc.before)


cat(
  "QC前细胞数：",
  ncol(gse292515),
  "\n"
)

gse292515 <- subset(
  gse292515,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

cat(
  "QC后细胞数：",
  ncol(gse292515),
  "\n"
)




p.qc.after <- VlnPlot(
  gse292515,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  pt.size = 0,
  ncol = 3,
  layer = "counts"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.title.x = element_blank()
  )

print(p.qc.after)




##############################################################
## 9. 保存QC后的对象
##############################################################

output.dir <- "results_Ism1_epithelial"

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  gse292515,
  file = file.path(
    output.dir,
    "GSE292515_after_QC.rds"
  )
)


##############################################################
## 10. 合并Seurat v5的RNA layers
##############################################################

DefaultAssay(gse292515) <- "RNA"

Layers(gse292515[["RNA"]])

gse292515 <- JoinLayers(
  gse292515,
  assay = "RNA"
)

Layers(gse292515[["RNA"]])








##############################################################
## 11. Normalize
##############################################################

gse292515 <- NormalizeData(
  gse292515,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


##############################################################
## 12. 找高变基因
##############################################################

gse292515 <- FindVariableFeatures(
  gse292515,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

length(
  VariableFeatures(gse292515)
)


##############################################################
## 13. Scale高变基因
##############################################################

gse292515 <- ScaleData(
  gse292515,
  features = VariableFeatures(gse292515),
  verbose = TRUE
)


##############################################################
## 14. PCA
##############################################################

gse292515 <- RunPCA(
  gse292515,
  features = VariableFeatures(gse292515),
  npcs = 30,
  verbose = TRUE
)


##############################################################
## 15. ElbowPlot
##############################################################

p.elbow <- ElbowPlot(
  gse292515,
  ndims = 30
)

print(p.elbow)

ggsave(
  filename = file.path(
    output.dir,
    "01_ElbowPlot.png"
  ),
  plot = p.elbow,
  width = 7,
  height = 5,
  dpi = 300
)



##############################################################
## 16. 设定后续使用的维度
##############################################################

dims.use <- 1:20





##############################################################
## 17. Harmony
##############################################################

library(harmony)

gse292515 <- harmony::RunHarmony(
  object = gse292515,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = dims.use,
  reduction.save = "harmony"
)

Reductions(gse292515)



##############################################################
## 18. Neighbor
##############################################################

gse292515 <- FindNeighbors(
  gse292515,
  reduction = "harmony",
  dims = dims.use
)


##############################################################
## 19. Cluster
##############################################################

gse292515 <- FindClusters(
  gse292515,
  resolution = 0.5
)


##############################################################
## 20. UMAP
##############################################################

gse292515 <- RunUMAP(
  gse292515,
  reduction = "harmony",
  dims = dims.use,
  umap.method = "uwot",
  metric = "cosine",
  seed.use = 1234
)

##############################################################
## 21. t-SNE
##############################################################

gse292515 <- RunTSNE(
  gse292515,
  reduction = "harmony",
  dims = dims.use,
  seed.use = 1234
)






##############################################################
## 22. Cluster UMAP
##############################################################

p.cluster <- DimPlot(
  object = gse292515,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.25,
  raster = TRUE,
  raster.dpi = c(1024, 1024)
) +
  ggtitle("GSE292515 clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

print(p.cluster)

ggsave(
  filename = file.path(
    output.dir,
    "02_UMAP_clusters.png"
  ),
  plot = p.cluster,
  width = 11,
  height = 8,
  dpi = 300
)




table(
  gse292515$seurat_clusters
)

length(
  unique(gse292515$seurat_clusters)
)


##############################################################
## 23. 按实验组画UMAP
##############################################################

gse292515$group <- factor(
  gse292515$group,
  levels = c(
    "Young_Mock",
    "Young_IAV",
    "Old_Mock",
    "Old_IAV"
  )
)

p.group <- DimPlot(
  object = gse292515,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.25,
  shuffle = TRUE,
  seed = 1234,
  raster = TRUE,
  raster.dpi = c(1024, 1024)
) +
  ggtitle("GSE292515 experimental groups") +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

print(p.group)

ggsave(
  filename = file.path(
    output.dir,
    "03_UMAP_groups.png"
  ),
  plot = p.group,
  width = 11,
  height = 8,
  dpi = 300
)


##############################################################
## 24. 按样本画UMAP
##############################################################

p.sample <- DimPlot(
  object = gse292515,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.25,
  shuffle = TRUE,
  seed = 1234,
  raster = TRUE,
  raster.dpi = c(1024, 1024)
) +
  ggtitle("GSE292515 samples") +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

print(p.sample)




##############################################################
## 25. 分组UMAP
##############################################################

p.group.split <- DimPlot(
  object = gse292515,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "group",
  label = FALSE,
  pt.size = 0.2,
  raster = TRUE,
  raster.dpi = c(800, 800),
  ncol = 2
) +
  NoLegend()

print(p.group.split)

ggsave(
  filename = file.path(
    output.dir,
    "04_UMAP_clusters_split_by_group.png"
  ),
  plot = p.group.split,
  width = 15,
  height = 12,
  dpi = 300
)






##############################################################
## 26. 经典肺细胞marker
##############################################################

Idents(gse292515) <- "seurat_clusters"

lung.marker.list <- list(
  
  "Pan-epithelial" = c(
    "Epcam",
    "Krt8",
    "Krt18",
    "Krt19"
  ),
  
  "AT1" = c(
    "Ager",
    "Pdpn",
    "Hopx",
    "Cav1",
    "Rtkn2"
  ),
  
  "AT2" = c(
    "Sftpc",
    "Sftpb",
    "Sftpa1",
    "Abca3",
    "Slc34a2"
  ),
  
  "Club" = c(
    "Scgb1a1",
    "Scgb3a1",
    "Cyp2f2"
  ),
  
  "Ciliated" = c(
    "Foxj1",
    "Pifo",
    "Tppp3"
  ),
  
  "Fibroblast" = c(
    "Col1a1",
    "Col1a2",
    "Dcn",
    "Lum",
    "Pdgfra"
  ),
  
  "Endothelial" = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Egfl7"
  ),
  
  "Macrophage" = c(
    "Lyz2",
    "Adgre1",
    "Csf1r",
    "Ctss"
  ),
  
  "Alveolar macrophage" = c(
    "Marco",
    "Siglecf",
    "Fabp4",
    "Pparg"
  ),
  
  "Monocyte" = c(
    "Ly6c2",
    "Ccr2",
    "Fcgr3",
    "Lyz2"
  ),
  
  "Neutrophil" = c(
    "S100a8",
    "S100a9",
    "Ly6g",
    "Mmp8"
  ),
  
  "Dendritic cell" = c(
    "Flt3",
    "Itgax",
    "Clec10a",
    "H2-Ab1"
  ),
  
  "T cell" = c(
    "Cd3d",
    "Cd3e",
    "Trbc1"
  ),
  
  "NK cell" = c(
    "Nkg7",
    "Klrk1",
    "Klrb1c"
  ),
  
  "B cell" = c(
    "Cd79a",
    "Ms4a1",
    "Cd74"
  )
)





##############################################################
## 27. 经典marker DotPlot
##############################################################

all.marker.genes <- unlist(
  lung.marker.list,
  use.names = FALSE
)

duplicated.marker.genes <- unique(
  all.marker.genes[
    duplicated(all.marker.genes)
  ]
)

duplicated.marker.genes

"Ly6c2" %in% rownames(gse292515)

grep(
  "Ly6c",
  rownames(gse292515),
  value = TRUE,
  ignore.case = TRUE
)






p.marker <- DotPlot(
  object = gse292515,
  features = lung.marker.list,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Canonical lung cell markers") +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )
  )

print(p.marker)

ggsave(
  filename = file.path(
    output.dir,
    "05_Canonical_lung_markers_by_cluster.png"
  ),
  plot = p.marker,
  width = 20,
  height = 10,
  dpi = 300
)
##############################################################
## 只保留矩阵中存在的marker
##############################################################

lung.marker.list <- lapply(
  lung.marker.list,
  function(x) {
    intersect(
      x,
      rownames(gse292515)
    )
  }
)

lung.marker.list <- lung.marker.list[
  lengths(lung.marker.list) > 0
]

print(lung.marker.list)






##############################################################
## 28. 主要谱系FeaturePlot
##############################################################

major.marker.genes <- intersect(
  c(
    "Ptprc",    # immune
    "Epcam",    # epithelial
    "Sftpc",    # AT2
    "Ager",     # AT1
    "Scgb1a1",  # Club
    "Foxj1",    # Ciliated
    "Pecam1",   # Endothelial
    "Col1a1",   # Fibroblast
    "Lyz2"      # Myeloid
  ),
  rownames(gse292515)
)

install.packages("ggrastr")
library(ggrastr)


p.major.marker <- FeaturePlot(
  object = gse292515,
  features = major.marker.genes,
  reduction = "umap",
  order = TRUE,
  ncol = 3,
  raster = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99"
)

print(p.major.marker)

ggsave(
  filename = file.path(
    output.dir,
    "06_Major_lineage_marker_FeaturePlots.png"
  ),
  plot = p.major.marker,
  width = 15,
  height = 14,
  dpi = 300
)



##############################################################
## 29. cluster UMAP
##############################################################

p.cluster <- DimPlot(
  object = gse292515,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.4,
  raster = FALSE
) +
  ggtitle("GSE292515 clusters") +
  theme_classic()

print(p.cluster)

ggsave(
  filename = file.path(
    output.dir,
    "07_UMAP_clusters_labeled.png"
  ),
  plot = p.cluster,
  width = 10,
  height = 8,
  dpi = 300
)




##############################################################
## 29. 检查Ism1
##############################################################

"Ism1" %in% rownames(gse292515)

grep(
  "Ism1",
  rownames(gse292515),
  value = TRUE,
  ignore.case = TRUE
)


##############################################################
## 30. 全部细胞中的Ism1 FeaturePlot
##############################################################

p.Ism1.all <- FeaturePlot(
  object = gse292515,
  features = "Ism1",
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  raster = TRUE
) +
  ggtitle("Ism1 expression in mouse lung")

print(p.Ism1.all)


ggsave(
  filename = file.path(
    output.dir,
    "07_Ism1_FeaturePlot_all_cells.png"
  ),
  plot = p.Ism1.all,
  width = 9,
  height = 7,
  dpi = 300
)


##############################################################
## 31. 按实验组显示Ism1
##############################################################

p.Ism1.group <- FeaturePlot(
  object = gse292515,
  features = "Ism1",
  reduction = "umap",
  split.by = "group",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  raster = TRUE,
  ncol = 2
)

print(p.Ism1.group)

ggsave(
  filename = file.path(
    output.dir,
    "08_Ism1_FeaturePlot_by_group.png"
  ),
  plot = p.Ism1.group,
  width = 14,
  height = 11,
  dpi = 300
)




##############################################################
## 32. Ism1按cluster表达
##############################################################

p.Ism1.cluster.dot <- DotPlot(
  object = gse292515,
  features = "Ism1",
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Ism1 expression by cluster")

print(p.Ism1.cluster.dot)

p.Ism1.cluster.vln <- VlnPlot(
  object = gse292515,
  features = "Ism1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  RotatedAxis() +
  ggtitle("Ism1 expression by cluster")

print(p.Ism1.cluster.vln)



##############################################################
## 33. Ism1 cluster表达汇总
##############################################################

Ism1.cluster.data <- FetchData(
  gse292515,
  vars = c(
    "Ism1",
    "seurat_clusters",
    "sample",
    "age",
    "treatment",
    "group"
  )
)

Ism1.cluster.summary <- Ism1.cluster.data |>
  group_by(
    seurat_clusters
  ) |>
  summarise(
    n_cells = n(),
    
    mean_Ism1 = mean(
      Ism1,
      na.rm = TRUE
    ),
    
    median_Ism1 = median(
      Ism1,
      na.rm = TRUE
    ),
    
    pct_Ism1_positive = mean(
      Ism1 > 0,
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  ) |>
  arrange(
    desc(mean_Ism1)
  )

print(
  Ism1.cluster.summary,
  n = 100
)

write.csv(
  Ism1.cluster.summary,
  file = file.path(
    output.dir,
    "Ism1_expression_by_cluster.csv"
  ),
  row.names = FALSE
)



##############################################################
## 34. 提取Young小鼠
##############################################################

gse292515.young <- subset(
  gse292515,
  subset = age == "Young"
)

table(
  gse292515.young$sample,
  gse292515.young$treatment
)

p.Ism1.young <- FeaturePlot(
  object = gse292515.young,
  features = "Ism1",
  reduction = "umap",
  split.by = "treatment",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  raster = TRUE
)

print(p.Ism1.young)

ggsave(
  filename = file.path(
    output.dir,
    "09_Ism1_Young_Mock_vs_IAV.png"
  ),
  plot = p.Ism1.young,
  width = 13,
  height = 6,
  dpi = 300
)




##############################################################
## 35. 提取Old小鼠
##############################################################

gse292515.old <- subset(
  gse292515,
  subset = age == "Old"
)

## 检查Old小鼠各样本和处理组
table(
  gse292515.old$sample,
  gse292515.old$treatment
)


##############################################################
## Old小鼠：Mock vs IAV的Ism1 FeaturePlot
##############################################################

p.Ism1.old <- FeaturePlot(
  object = gse292515.old,
  features = "Ism1",
  reduction = "umap",
  split.by = "treatment",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  raster = TRUE
)

print(p.Ism1.old)


##############################################################
## 保存图片
##############################################################

ggsave(
  filename = file.path(
    output.dir,
    "10_Ism1_Old_Mock_vs_IAV.png"
  ),
  plot = p.Ism1.old,
  width = 13,
  height = 6,
  dpi = 300,
  bg = "white"
)



##############################################################
## Young/Old × Mock/IAV：四组统一展示
##############################################################

gse292515$group <- factor(
  gse292515$group,
  levels = c(
    "Young_Mock",
    "Young_IAV",
    "Old_Mock",
    "Old_IAV"
  )
)

p.Ism1.four.groups <- FeaturePlot(
  object = gse292515,
  features = "Ism1",
  reduction = "umap",
  split.by = "group",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  keep.scale = "all",
  raster = TRUE,
  ncol = 2
)

print(p.Ism1.four.groups)

ggsave(
  filename = file.path(
    output.dir,
    "11_Ism1_Young_Old_Mock_IAV.png"
  ),
  plot = p.Ism1.four.groups,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)




