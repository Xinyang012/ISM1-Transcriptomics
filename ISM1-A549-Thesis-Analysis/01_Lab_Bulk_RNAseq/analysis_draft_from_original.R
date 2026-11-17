# Thesis source-code archive — reviewer draft, NOT validated end-to-end.
# Original analysis logic, thresholds and sample identifiers preserved.
# Replace PATH_TO_* placeholders with permitted local paths before running.
# Check README.md and notes/VALIDATION_CHECKLIST.md before sharing publicly.

### 0. 安装 & 加载包
############################################################
install.packages(c("dplyr","tibble","ggplot2","pheatmap","msigdbr"))
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db","enrichplot"))

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

############################################################
### 1. 读取数据
############################################################
setwd("PATH_TO_LAB_RNASEQ_INPUTS")

df_3h <- read.csv("RNA_data_3h.csv", row.names = 1)
df_12h <- read.csv("RNA_data_12h.csv", row.names = 1)

# 去掉 ENSEMBL 版本号
rownames(df_3h) <- sub("\\..*","",rownames(df_3h))
rownames(df_12h) <- sub("\\..*","",rownames(df_12h))

############################################################
### 2. 样本分组信息
############################################################
samples_list <- list(
  "3h_SEN" = c("X3h_SEN1","X3h_SEN2","X3h_SEN3","X3h_SEN.ISM1_1","X3h_SEN.ISM1_2","X3h_SEN.ISM1_3"),
  "3h_NOR" = c("X3h_NOR1","X3h_NOR2","X3h_NOR3","X3h_NOR.ISM1_1","X3h_NOR.ISM1_2","X3h_NOR.ISM1_3"),
  "12h_SEN" = c("X12h_SEN1","X12h_SEN2","X12h_SEN3","X12h_SEN.ISM1_1","X12h_SEN.ISM1_2","X12h_SEN.ISM1_3"),
  "12h_NOR" = c("X12h_NOR1","X12h_NOR2","X12h_NOR3","X12h_NOR.ISM1_1","X12h_NOR.ISM1_2","X12h_NOR.ISM1_3")
)
df_list <- list("3h_SEN"=df_3h,"3h_NOR"=df_3h,"12h_SEN"=df_12h,"12h_NOR"=df_12h)

############################################################
### 3. MSigDB Hallmark gene sets
############################################################
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(x = hallmark_sets$entrez_gene, f = hallmark_sets$gs_name)


hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)  # TERM2GENE 格式
hallmark_list <- hallmark

############################################################
### 3. 批量 GSEA Hallmark 函数（修正版）
############################################################

run_hallmark_gsea <- function(name, sel_samples, df, hallmark_list=hallmark){
  
  message("Processing: ", name)
  
  # 分组信息
  group <- ifelse(grepl("ISM", sel_samples, ignore.case = TRUE),
                  paste0(name, ".ISM"), name)
  sample_info <- data.frame(
    sample = sel_samples,
    group = factor(group)
  )
  rownames(sample_info) <- sel_samples
  
  # count matrix
  count_mat <- df[, sel_samples, drop=FALSE]
  count_mat <- count_mat[rowSums(count_mat) > 0, ]
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = sample_info,
    design = ~ group
  )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", paste0(name, ".ISM"), name))
  res <- as.data.frame(res)
  res$ENSEMBL <- rownames(res)
  res <- res[!is.na(res$stat), ]
  
  # ENSEMBL → ENTREZID
  gene_map <- bitr(res$ENSEMBL, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]
  
  # 合并并去重，每个 ENTREZID 取 stat 最大值
  res2 <- merge(res, gene_map, by="ENSEMBL")
  res2 <- res2 %>%
    group_by(ENTREZID) %>%
    summarise(stat = max(stat)) %>%
    ungroup()
  
  # 构建安全 numeric 向量
  geneList <- res2$stat
  names(geneList) <- res2$ENTREZID
  geneList <- sort(as.numeric(geneList), decreasing=TRUE)
  names(geneList) <- res2$ENTREZID
  
  # GSEA Hallmark
  gsea_res <- GSEA(geneList, TERM2GENE=hallmark_list, pvalueCutoff=0.05)
  
  # 提取结果并画 Top10 条形图
  gsea_top <- gsea_res@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
    group_by(Direction) %>%
    slice_max(order_by = abs(NES), n=10) %>%
    ungroup()
  
  gsea_top$Description <- factor(gsea_top$Description,
                                 levels = rev(unique(gsea_top$Description)))
  
  p <- ggplot(gsea_top, aes(x = NES, y = Description, fill = Direction)) +
    geom_col() +
    facet_wrap(~Direction, scales = "free_y") +
    scale_fill_manual(values = c("Upregulated"="red","Downregulated"="blue")) +
    theme_bw() +
    labs(x="Normalized Enrichment Score (NES)",
         y="Pathway",
         title=paste0("GSEA Hallmark NES: ", name)) +
    theme(strip.text = element_text(size=12, face="bold"),
          axis.text.y = element_text(size=10),
          legend.position="none")
  
  # 保存图片
  ggsave(paste0("GSEA_Hallmark_", name, ".png"), p, width=16, height=6, dpi=300)
  
  return(gsea_res)
}

############################################################
### 4. 批量运行
############################################################

gsea_results <- list()
for(name in names(samples_list)){
  gsea_results[[name]] <- run_hallmark_gsea(name, samples_list[[name]], df_list[[name]])
}

############################################################
### 完成 ✅





##########
#我现在又重新跑一下看看吧

library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

run_hallmark_gsea_fixed <- function(name, samples, df){
  
  cat("Processing:", name, "\n")
  
  # -----------------------------
  # 分组信息
  # -----------------------------
  group <- ifelse(grepl("ISM", samples, ignore.case = TRUE),
                  paste0(sub("_.*", "", name), ".ISM1"),
                  sub("_.*", "", name))
  
  sample_info <- data.frame(sample = samples, group = factor(group))
  rownames(sample_info) <- samples
  
  # -----------------------------
  # count matrix
  # -----------------------------
  count_mat <- df[, samples]
  count_mat <- count_mat[rowSums(count_mat) > 0, ]
  
  # -----------------------------
  # DESeq2
  # -----------------------------
  dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                colData = sample_info,
                                design = ~ group)
  
  dds <- DESeq(dds)
  res <- results(dds,
                 contrast = c("group", paste0(sub("_.*", "", name), ".ISM1"),
                              sub("_.*", "", name)))
  res <- as.data.frame(res)
  res$ENSEMBL <- rownames(res)
  res <- res[!is.na(res$stat), ]
  
  # -----------------------------
  # geneList numeric
  # -----------------------------
  gene_list <- res$stat
  names(gene_list) <- res$ENSEMBL
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # -----------------------------
  # ENSEMBL -> ENTREZID
  # -----------------------------
  gene_map <- bitr(names(gene_list), fromType = "ENSEMBL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # 去重复
  gene_list2 <- gene_list[gene_map$ENSEMBL]
  df2 <- data.frame(ENTREZID = gene_map$ENTREZID, stat = gene_list2)
  df2 <- df2 %>%
    group_by(ENTREZID) %>%
    summarise(stat = max(stat))
  
  gene_list_final <- as.numeric(df2$stat)
  names(gene_list_final) <- df2$ENTREZID
  gene_list_final <- sort(gene_list_final, decreasing = TRUE)
  
  # -----------------------------
  # GSEA Hallmark
  # -----------------------------
  gsea_res <- GSEA(gene_list_final,
                   TERM2GENE = msigdbr(species = "Homo sapiens", category = "H"),
                   pvalueCutoff = 0.05,
                   verbose = FALSE)
  
  # -----------------------------
  # 提取结果 & 上下调
  # -----------------------------
  gsea_df <- gsea_res@result
  gsea_df <- gsea_df %>%
    filter(p.adjust < 0.05) %>%
    mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
    mutate(Description = gsub("^HALLMARK_", "", Description))
  
  # Top10 上下调
  gsea_top <- gsea_df %>%
    group_by(Direction) %>%
    slice_max(order_by = abs(NES), n = 10)
  gsea_top$Description <- factor(gsea_top$Description,
                                 levels = rev(unique(gsea_top$Description)))
  
  # -----------------------------
  # 画条形图
  # -----------------------------
  p <- ggplot(gsea_top, aes(x = NES, y = Description, fill = Direction)) +
    geom_col() +
    facet_wrap(~Direction, scales = "free_y") +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    theme_bw() +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway",
         title = paste0("GSEA Hallmark: ", name)) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  
  # 保存图
  ggsave(paste0("GSEA_Hallmark_", name, ".png"),
         plot = p,
         width = 14,
         height = 10,   # 高度加大避免纵坐标挤压
         dpi = 300)
  
  return(gsea_res)
}

library(msigdbr)  # 确保加载

gsea_results <- list()
for(name in names(samples_list)){
  gsea_results[[name]] <- run_hallmark_gsea_fixed(name, samples_list[[name]], df_list[[name]])
}



#####
#我又要重新跑了
############################################################
# 设置工作目录
############################################################
setwd("PATH_TO_LAB_RNASEQ_INPUTS")

############################################################
# 设置工作目录
############################################################
setwd("PATH_TO_LAB_RNASEQ_INPUTS")

############################################################
# 读取数据
############################################################
df_3h <- read.csv("RNA_data_3h.csv", row.names = 1)
df_12h <- read.csv("RNA_data_12h.csv", row.names = 1)

# 去掉 ENSEMBL 版本号
rownames(df_3h) <- sub("\\..*","",rownames(df_3h))
rownames(df_12h) <- sub("\\..*","",rownames(df_12h))

############################################################
# 样本分组信息
############################################################
samples_list <- list(
  "3h_SEN" = c("X3h_SEN1","X3h_SEN2","X3h_SEN3","X3h_SEN.ISM1_1","X3h_SEN.ISM1_2","X3h_SEN.ISM1_3"),
  "3h_NOR" = c("X3h_NOR1","X3h_NOR2","X3h_NOR3","X3h_NOR.ISM1_1","X3h_NOR.ISM1_2","X3h_NOR.ISM1_3"),
  "12h_SEN" = c("X12h_SEN1","X12h_SEN2","X12h_SEN3","X12h_SEN.ISM1_1","X12h_SEN.ISM1_2","X12h_SEN.ISM1_3"),
  "12h_NOR" = c("X12h_NOR1","X12h_NOR2","X12h_NOR3","X12h_NOR.ISM1_1","X12h_NOR.ISM1_2","X12h_NOR.ISM1_3")
)

df_list <- list("3h_SEN"=df_3h,"3h_NOR"=df_3h,"12h_SEN"=df_12h,"12h_NOR"=df_12h)

############################################################
# 安装必要包（如未安装）
############################################################
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db","msigdbr"))

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

############################################################
# 获取 Hallmark gene sets
############################################################
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

############################################################
samples_list <- list(
  "3h_SEN" = c("X3h_SEN1","X3h_SEN2","X3h_SEN3"),
  "3h_SEN_ISM1" = c("X3h_SEN.ISM1_1","X3h_SEN.ISM1_2","X3h_SEN.ISM1_3"),
  
  "3h_NOR" = c("X3h_NOR1","X3h_NOR2","X3h_NOR3"),
  "3h_NOR_ISM1" = c("X3h_NOR.ISM1_1","X3h_NOR.ISM1_2","X3h_NOR.ISM1_3"),
  
  "12h_SEN" = c("X12h_SEN1","X12h_SEN2","X12h_SEN3"),
  "12h_SEN_ISM1" = c("X12h_SEN.ISM1_1","X12h_SEN.ISM1_2","X12h_SEN.ISM1_3"),
  
  "12h_NOR" = c("X12h_NOR1","X12h_NOR2","X12h_NOR3"),
  "12h_NOR_ISM1" = c("X12h_NOR.ISM1_1","X12h_NOR.ISM1_2","X12h_NOR.ISM1_3")
)
comparisons <- list(
  c("3h_SEN_ISM1","3h_SEN"),
  c("3h_NOR_ISM1","3h_NOR"),
  c("12h_SEN_ISM1","12h_SEN"),
  c("12h_NOR_ISM1","12h_NOR")
)
gsea_results <- list()



run_DESeq2 <- function(counts, group1, group2){
  
  group1 <- intersect(group1, colnames(counts))
  group2 <- intersect(group2, colnames(counts))
  
  if(length(group1)==0 | length(group2)==0){
    stop("样本名不匹配，请检查！")
  }
  
  sub_counts <- counts[, c(group1, group2)]
  
  colData <- data.frame(
    row.names = colnames(sub_counts),
    condition = factor(c(rep("ISM1", length(group1)),
                         rep("CTRL", length(group2))))
  )
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts,
                                colData = colData,
                                design = ~ condition)
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  res <- res[order(res$padj), ]
  res <- res[!is.na(res$log2FoldChange), ]
  
  gene_list <- res$log2FoldChange
  names(gene_list) <- rownames(res)
  
  # 打破重复值（避免GSEA警告）
  gene_list <- gene_list + rnorm(length(gene_list), 0, 1e-6)
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}


convert_to_symbol <- function(gene_list){
  
  gene_df <- bitr(names(gene_list),
                  fromType = "ENSEMBL",
                  toType = "SYMBOL",
                  OrgDb = org.Hs.eg.db)
  
  gene_list <- gene_list[gene_df$ENSEMBL]
  names(gene_list) <- gene_df$SYMBOL
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}


plot_gsea_bar <- function(gsea_obj, top_n = 10, title = ""){
  
  df <- as.data.frame(gsea_obj@result)
  
  if(nrow(df)==0){
    return(NULL)
  }
  
  df <- df %>% filter(!is.na(NES))
  
  df_up <- df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_n)
  df_down <- df %>% filter(NES < 0) %>% arrange(NES) %>% head(top_n)
  
  df_up$Direction <- "Upregulated"
  df_down$Direction <- "Downregulated"
  
  plot_df <- rbind(df_up, df_down)
  
  plot_df$NES_abs <- abs(plot_df$NES)
  
  plot_df$Description <- gsub("HALLMARK_", "", plot_df$Description)
  plot_df$Description <- gsub("_", " ", plot_df$Description)
  
  plot_df$Description <- factor(plot_df$Description,
                                levels = rev(unique(plot_df$Description)))
  
  ggplot(plot_df, aes(x = NES_abs, y = Description, fill = Direction)) +
    geom_bar(stat="identity") +
    facet_wrap(~Direction, scales="free_y") +
    scale_fill_manual(values=c("Upregulated"="red","Downregulated"="blue")) +
    labs(title=title, x="|NES|", y="Pathway") +
    theme_bw() +
    theme(
      strip.text = element_text(angle=90, size=12),
      plot.title = element_text(size=14, face="bold"),
      legend.position="none"
    )
}


comparisons <- list(
  c("3h_SEN_ISM1","3h_SEN"),
  c("3h_NOR_ISM1","3h_NOR"),
  c("12h_SEN_ISM1","12h_SEN"),
  c("12h_NOR_ISM1","12h_NOR")
)

gsea_results <- list()

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_df <- hallmark[, c("gs_name", "gene_symbol")]

run_DESeq2 <- function(counts, group1, group2){
  
  group1 <- intersect(group1, colnames(counts))
  group2 <- intersect(group2, colnames(counts))
  
  sub_counts <- counts[, c(group1, group2)]
  
  colData <- data.frame(
    row.names = colnames(sub_counts),
    condition = factor(c(rep("ISM1", length(group1)),
                         rep("CTRL", length(group2))))
  )
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts,
                                colData = colData,
                                design = ~ condition)
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  res <- res[!is.na(res$log2FoldChange), ]
  
  gene_list <- res$log2FoldChange
  names(gene_list) <- rownames(res)
  
  gene_list <- gene_list + rnorm(length(gene_list), 0, 1e-6)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

convert_to_symbol <- function(gene_list){
  
  gene_df <- bitr(names(gene_list),
                  fromType = "ENSEMBL",
                  toType = "SYMBOL",
                  OrgDb = org.Hs.eg.db)
  
  gene_df$logFC <- gene_list[gene_df$ENSEMBL]
  
  gene_df <- gene_df[!is.na(gene_df$SYMBOL), ]
  
  # ⭐按变化幅度排序 + 去重
  gene_df <- gene_df[order(abs(gene_df$logFC), decreasing = TRUE), ]
  gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]
  
  gene_list2 <- gene_df$logFC
  names(gene_list2) <- gene_df$SYMBOL
  
  gene_list2 <- sort(gene_list2, decreasing = TRUE)
  
  return(gene_list2)
}

plot_gsea_bar <- function(gsea_obj, title=""){
  
  df <- as.data.frame(gsea_obj@result)
  
  if(nrow(df)==0){
    return(NULL)
  }
  
  df_up <- df[df$NES > 0, ]
  df_down <- df[df$NES < 0, ]
  
  df_up <- df_up[order(-df_up$NES), ][1:min(10,nrow(df_up)), ]
  df_down <- df_down[order(df_down$NES), ][1:min(10,nrow(df_down)), ]
  
  df_up$Direction <- "Up"
  df_down$Direction <- "Down"
  
  plot_df <- rbind(df_up, df_down)
  
  plot_df$NES <- abs(plot_df$NES)
  plot_df$Description <- gsub("HALLMARK_", "", plot_df$Description)
  
  ggplot(plot_df, aes(NES, reorder(Description, NES), fill=Direction)) +
    geom_col() +
    facet_wrap(~Direction, scales="free_y") +
    scale_fill_manual(values=c("Up"="red","Down"="blue")) +
    theme_bw() +
    labs(title=title, x="|NES|", y="Pathway")
}


run_one <- function(comp){
  
  group1 <- samples_list[[comp[1]]]
  group2 <- samples_list[[comp[2]]]
  
  df_use <- if(grepl("3h", comp[1])) df_3h else df_12h
  
  gene_list <- run_DESeq2(df_use, group1, group2)
  gene_list <- convert_to_symbol(gene_list)
  
  gsea <- tryCatch({
    GSEA(geneList = gene_list,
         TERM2GENE = hallmark_df,
         pvalueCutoff = 0.05)
  }, error=function(e){
    cat("GSEA error:", comp, "\n")
    return(NULL)
  })
  
  if(is.null(gsea)) return(NULL)
  
  p <- plot_gsea_bar(gsea, paste(comp, collapse="_vs_"))
  
  if(!is.null(p)){
    print(p)
    ggsave(paste0(paste(comp, collapse="_vs_"), ".png"),
           p, width=6, height=8)
  } else {
    cat("No enrichment:", comp, "\n")
  }
}

run_one(c("3h_SEN_ISM1","3h_SEN"))
run_one(c("3h_NOR_ISM1","3h_NOR"))
run_one(c("12h_SEN_ISM1","12h_SEN"))
run_one(c("12h_NOR_ISM1","12h_NOR"))