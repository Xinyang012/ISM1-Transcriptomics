# Recovered source: prior user-uploaded new.R; provenance retained.
# NOTE: requires objects created before this script; see README.
# NOT independently executed/verified.

############################################################
## Figure 3B
## Ism1-associated genes in IAV-infected lungs
############################################################

library(DESeq2)
library(dplyr)
library(ggplot2)

############################################################
## 1. 检查现有对象
############################################################

stopifnot(
  exists("vst.matrix.iav271812"),
  exists("metadata.timecourse271812")
)

cat("VST matrix dimensions:\n")
print(dim(vst.matrix.iav271812))

cat("\nIAV metadata:\n")

meta.iav271812 <- metadata.timecourse271812 |>
  dplyr::filter(
    infection_status == "Influenza"
  )

print(
  meta.iav271812[
    ,
    c(
      "sample",
      "time_group",
      "day",
      "replicate"
    )
  ],
  row.names = FALSE
)

stopifnot(
  nrow(meta.iav271812) == 12
)

############################################################
## 2. 确认 sample 顺序
############################################################

meta.iav271812 <- meta.iav271812[
  match(
    colnames(vst.matrix.iav271812),
    meta.iav271812$sample
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(vst.matrix.iav271812),
    meta.iav271812$sample
  )
)

############################################################
## 3. Ism1 expression
############################################################

stopifnot(
  "Ism1" %in%
    rownames(vst.matrix.iav271812)
)

ism1.expression.iav271812 <-
  as.numeric(
    vst.matrix.iav271812[
      "Ism1",
    ]
  )

############################################################
## 4. Spearman correlation
############################################################

cor.results.iav271812 <- lapply(
  
  setdiff(
    rownames(vst.matrix.iav271812),
    "Ism1"
  ),
  
  function(gene) {
    
    gene.expression <-
      as.numeric(
        vst.matrix.iav271812[
          gene,
        ]
      )
    
    test <-
      suppressWarnings(
        cor.test(
          ism1.expression.iav271812,
          gene.expression,
          method = "spearman",
          exact = FALSE
        )
      )
    
    data.frame(
      
      gene = gene,
      
      rho =
        unname(
          test$estimate
        ),
      
      pvalue =
        test$p.value,
      
      stringsAsFactors = FALSE
    )
  }
)

############################################################
## 5. Combine
############################################################

cor.results.iav271812 <-
  dplyr::bind_rows(
    cor.results.iav271812
  )

############################################################
## 6. BH correction
############################################################

cor.results.iav271812$padj <-
  p.adjust(
    cor.results.iav271812$pvalue,
    method = "BH"
  )

############################################################
## 7. Rank
############################################################

cor.results.iav271812 <-
  cor.results.iav271812 |>
  dplyr::arrange(
    desc(rho)
  )

############################################################
## 8. 保存
############################################################

write.csv(
  cor.results.iav271812,
  file = file.path(
    result.dir.271812,
    "GSE271812_Ism1_IAV_Spearman_correlation.csv"
  ),
  row.names = FALSE
)

############################################################
## 9. 查看结果
############################################################

cat("\n==============================\n")
cat("Top positive Ism1-associated genes\n")
cat("==============================\n")

print(
  cor.results.iav271812 |>
    dplyr::filter(
      rho > 0
    ) |>
    dplyr::arrange(
      padj
    ) |>
    head(30),
  row.names = FALSE
)

cat("\n==============================\n")
cat("Top negative Ism1-associated genes\n")
cat("==============================\n")

print(
  cor.results.iav271812 |>
    dplyr::filter(
      rho < 0
    ) |>
    dplyr::arrange(
      padj
    ) |>
    head(30),
  row.names = FALSE
)
















############################################################
## Figure 3B
## Sensitivity analysis:
## association after accounting for infection day
############################################################

############################################################
## 1. 构建 IAV expression data
############################################################

day.iav271812 <-
  as.numeric(
    meta.iav271812$day
  )

############################################################
## 2. Ism1 residual after removing day effect
############################################################

ism1.residual271812 <-
  residuals(
    lm(
      ism1.expression.iav271812 ~
        day.iav271812
    )
  )

############################################################
## 3. Calculate day-adjusted correlations
############################################################

partial.cor.results.iav271812 <-
  lapply(
    
    setdiff(
      rownames(
        vst.matrix.iav271812
      ),
      "Ism1"
    ),
    
    function(gene) {
      
      gene.expression <-
        as.numeric(
          vst.matrix.iav271812[
            gene,
          ]
        )
      
      gene.residual <-
        residuals(
          lm(
            gene.expression ~
              day.iav271812
          )
        )
      
      test <-
        suppressWarnings(
          cor.test(
            ism1.residual271812,
            gene.residual,
            method = "spearman",
            exact = FALSE
          )
        )
      
      data.frame(
        
        gene = gene,
        
        rho_day_adjusted =
          unname(
            test$estimate
          ),
        
        pvalue_day_adjusted =
          test$p.value,
        
        stringsAsFactors = FALSE
      )
    }
  )

partial.cor.results.iav271812 <-
  dplyr::bind_rows(
    partial.cor.results.iav271812
  )

partial.cor.results.iav271812$padj_day_adjusted <-
  p.adjust(
    partial.cor.results.iav271812$pvalue_day_adjusted,
    method = "BH"
  )

partial.cor.results.iav271812 <-
  partial.cor.results.iav271812 |>
  dplyr::arrange(
    desc(rho_day_adjusted)
  )

############################################################
## 4. Save
############################################################

write.csv(
  partial.cor.results.iav271812,
  file = file.path(
    result.dir.271812,
    "GSE271812_Ism1_IAV_Spearman_day_adjusted.csv"
  ),
  row.names = FALSE
)

############################################################
## 5. 查看
############################################################

cat(
  "\nTop positive associations after adjusting for day:\n"
)

print(
  partial.cor.results.iav271812 |>
    dplyr::filter(
      rho_day_adjusted > 0
    ) |>
    dplyr::arrange(
      padj_day_adjusted
    ) |>
    head(30),
  row.names = FALSE
)




cat(
  "\nTop negative associations after adjusting for day:\n"
)

print(
  partial.cor.results.iav271812 |>
    dplyr::filter(
      rho_day_adjusted < 0
    ) |>
    dplyr::arrange(
      padj_day_adjusted
    ) |>
    head(30),
  row.names = FALSE
)





library(ggplot2)

plot.data <- data.frame(
  Ism1 = as.numeric(expr.iav271812["Ism1", ]),
  Ehd3 = as.numeric(expr.iav271812["Ehd3", ]),
  sample = colnames(expr.iav271812)
)

plot.data$day <- factor(
  sub(".*D(4|6|8).*", "\\1", plot.data$sample),
  levels = c("4", "6", "8")
)

plot.data


ls()
ls(pattern = "iav271812")
ls(pattern = "expr")


c(
  "Ism1" %in% rownames(expr.matrix.iav271812),
  "Ehd3" %in% rownames(expr.matrix.iav271812)
)




plot.data <- data.frame(
  Ism1 = as.numeric(expr.matrix.iav271812["Ism1", ]),
  Ehd3 = as.numeric(expr.matrix.iav271812["Ehd3", ]),
  sample = colnames(expr.matrix.iav271812)
)

plot.data




library(ggplot2)

plot.data$day <- day.iav271812

ggplot(
  plot.data,
  aes(
    x = Ism1,
    y = Ehd3,
    color = factor(day)
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  theme_classic() +
  labs(
    x = "Ism1 expression",
    y = "Ehd3 expression",
    color = "Day"
  )






plot.data$day <- day.iav271812

plot.data





day.iav271812
meta.iav271812[, c("geo_accession", "day")]
colnames(meta.iav271812)
colnames(expr.matrix.iav271812)



library(ggplot2)

ggplot(
  plot.data,
  aes(
    x = Ism1,
    y = Ehd3,
    color = factor(day)
  )
) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "Ism1 expression",
    y = "Ehd3 expression",
    color = "Day"
  )




cor.test(
  plot.data$Ism1,
  plot.data$Ehd3,
  method = "spearman",
  exact = FALSE
)









library(ggplot2)

ism1.plot.data <- data.frame(
  Ism1 = as.numeric(ism1.expression.iav271812),
  Day = factor(day.iav271812, levels = c(4, 6, 8))
)

ggplot(ism1.plot.data, aes(x = Day, y = Ism1)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.08,
    size = 3
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Day post infection",
    y = "ISM1 expression",
    title = "ISM1 expression during IAV infection"
  )





cor.test(
  plot.data$Ism1,
  plot.data$Ehd3,
  method = "spearman",
  exact = FALSE
)



plot.data <- data.frame(
  Ism1 = as.numeric(expr.matrix.iav271812["Ism1", ]),
  Ehd3 = as.numeric(expr.matrix.iav271812["Ehd3", ]),
  Day = day.iav271812,
  sample = colnames(expr.matrix.iav271812)
)


plot.data$Ism1_resid <- resid(
  lm(Ism1 ~ Day, data = plot.data)
)

plot.data$Ehd3_resid <- resid(
  lm(Ehd3 ~ Day, data = plot.data)
)




plot.data


ggplot(
  plot.data,
  aes(x = Ism1_resid, y = Ehd3_resid)
) +
  geom_point(
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "ISM1 expression (Day-adjusted residual)",
    y = "Ehd3 expression (Day-adjusted residual)",
    title = "Ehd3 is positively associated with ISM1 after adjustment for day"
  )


cor.test(
  plot.data$Ism1_resid,
  plot.data$Ehd3_resid,
  method = "pearson"
)




genes.of.interest <- c(
  "Ism1",
  "Ehd3",
  "Slfn14",
  "Gm41055"
)

candidate.plot.data <- data.frame(
  sample = colnames(expr.matrix.iav271812),
  Day = factor(day.iav271812, levels = c(4, 6, 8))
)

for (g in genes.of.interest) {
  candidate.plot.data[[g]] <-
    as.numeric(expr.matrix.iav271812[g, ])
}




library(tidyr)
library(dplyr)

candidate.long <- candidate.plot.data |>
  pivot_longer(
    cols = all_of(genes.of.interest),
    names_to = "Gene",
    values_to = "Expression"
  )


ggplot(
  candidate.long,
  aes(x = Day, y = Expression)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.5
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 2
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Day post infection",
    y = "Expression"
  )








make.partial.plot.data <- function(
    expr.matrix,
    gene,
    day
) {
  
  dat <- data.frame(
    Ism1 = as.numeric(expr.matrix["Ism1", ]),
    Gene = as.numeric(expr.matrix[gene, ]),
    Day = day,
    sample = colnames(expr.matrix)
  )
  
  dat$Ism1_resid <- resid(
    lm(Ism1 ~ Day, data = dat)
  )
  
  dat$Gene_resid <- resid(
    lm(Gene ~ Day, data = dat)
  )
  
  dat
}



ehd3.data <- make.partial.plot.data(
  expr.matrix.iav271812,
  "Ehd3",
  day.iav271812
)

slfn14.data <- make.partial.plot.data(
  expr.matrix.iav271812,
  "Slfn14",
  day.iav271812
)

gm41055.data <- make.partial.plot.data(
  expr.matrix.iav271812,
  "Gm41055",
  day.iav271812
)




ggplot(
  slfn14.data,
  aes(x = Ism1_resid, y = Gene_resid)
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "ISM1 expression\n(Day-adjusted)",
    y = "Slfn14 expression\n(Day-adjusted)",
    title = "Slfn14"
  )



ggplot(
  gm41055.data,
  aes(x = Ism1_resid, y = Gene_resid)
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "ISM1 expression\n(Day-adjusted)",
    y = "Gm41055 expression\n(Day-adjusted)",
    title = "Gm41055"
  )






library(pheatmap)

heat.genes <- c(
  "Ism1",
  "Ehd3",
  "Slfn14",
  "Gm41055"
)

heat.mat <- expr.matrix.iav271812[
  heat.genes,
  order(day.iav271812)
]

heat.mat.z <- t(
  scale(t(heat.mat))
)

annotation.col <- data.frame(
  Day = factor(
    day.iav271812[
      order(day.iav271812)
    ],
    levels = c(4, 6, 8)
  )
)

rownames(annotation.col) <-
  colnames(heat.mat.z)

pheatmap(
  heat.mat.z,
  annotation_col = annotation.col,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 12,
  main = "ISM1-associated genes during IAV infection"
)






head(partial.cor.results.iav271812)
dim(partial.cor.results.iav271812)
summary(partial.cor.results.iav271812$rho_day_adjusted)
sum(
  !is.na(partial.cor.results.iav271812$rho_day_adjusted)
)







n <- 12
k <- 1


partial.cor.results.iav271812$t_stat <- with(
  partial.cor.results.iav271812,
  rho_day_adjusted *
    sqrt((n - k - 2) / (1 - rho_day_adjusted^2))
)






gsea.rank.iav271812 <- partial.cor.results.iav271812 |>
  dplyr::filter(
    !is.na(gene),
    gene != "",
    !is.na(t_stat),
    is.finite(t_stat)
  ) |>
  dplyr::arrange(desc(t_stat))




head(gsea.rank.iav271812, 20)
tail(gsea.rank.iav271812, 20)






library(msigdbr)
library(clusterProfiler)
library(dplyr)





hallmark.mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)



hallmark.list <- split(
  hallmark.mouse$gene_symbol,
  hallmark.mouse$gs_name
)




gsea.rank <- gsea.rank.iav271812$t_stat
names(gsea.rank) <- gsea.rank.iav271812$gene




gsea.rank <- gsea.rank[
  !duplicated(names(gsea.rank))
]



gsea.rank <- sort(
  gsea.rank,
  decreasing = TRUE
)



library(clusterProfiler)

gsea.hallmark.ism1 <- GSEA(
  geneList = gsea.rank,
  TERM2GENE = hallmark.mouse[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)




as.data.frame(gsea.hallmark.ism1)




gsea.df <- as.data.frame(gsea.hallmark.ism1)


gsea.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(desc(NES)) |>
  head(15)



gsea.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(NES) |>
  head(15)


dot.df <- gsea.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(desc(abs(NES))) |>
  head(15)


library(ggplot2)

ggplot(
  dot.df,
  aes(
    x = NES,
    y = reorder(Description, NES),
    size = -log10(p.adjust),
    color = NES
  )
) +
  geom_point() +
  theme_classic(base_size = 12) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10 adjusted P",
    color = "NES",
    title = "Hallmark pathways associated with ISM1 expression"
  )






library(org.Mm.eg.db)
gsea.go.ism1 <- gseGO(
  geneList = gsea.rank,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)
gsea.go.df <- as.data.frame(gsea.go.ism1)
gsea.go.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(desc(NES)) |>
  head(20)
gsea.go.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(NES) |>
  head(20)








gsea.hallmark.df <- as.data.frame(gsea.hallmark.ism1)

gsea.hallmark.sig <- gsea.hallmark.df |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(desc(NES))

gsea.hallmark.sig



gsea.hallmark.sig |>
  filter(NES > 0) |>
  head(20)


gsea.hallmark.sig |>
  filter(NES < 0) |>
  arrange(NES) |>
  head(20)




gsea.go.df <- as.data.frame(gsea.go.ism1)

gsea.go.sig <- gsea.go.df |>
  filter(p.adjust < 0.05) |>
  arrange(desc(NES))

gsea.go.sig


gsea.go.sig |>
  filter(
    grepl(
      "JAK|STAT|coagulation|hemostasis|platelet|thromb|interferon|TLR|inflamm",
      Description,
      ignore.case = TRUE
    )
  )





library(dplyr)
library(ggplot2)

hallmark.df <- as.data.frame(gsea.hallmark.ism1)

hallmark.sig <- hallmark.df %>%
  filter(p.adjust < 0.05)

hallmark.up <- hallmark.sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(10)

hallmark.down <- hallmark.sig %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(10)

hallmark.plot.df <- bind_rows(
  hallmark.up,
  hallmark.down
) %>%
  arrange(NES)

ggplot(
  hallmark.plot.df,
  aes(
    x = NES,
    y = reorder(Description, NES),
    size = -log10(p.adjust)
  )
) +
  geom_point() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10 adjusted P"
  )



go.sig <- gsea.go.df %>%
  filter(p.adjust < 0.05)

go.up <- go.sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(15)

go.down <- go.sig %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(15)

go.plot.df <- bind_rows(
  go.up,
  go.down
) %>%
  arrange(NES)

ggplot(
  go.plot.df,
  aes(
    x = NES,
    y = reorder(Description, NES),
    size = -log10(p.adjust)
  )
) +
  geom_point() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10 adjusted P"
  )







library(dplyr)
library(ggplot2)
library(stringr)

kegg.sig <- kegg.df %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    direction = ifelse(NES > 0, "Positive", "Negative"),
    pathway = str_wrap(Description, width = 38)
  )

# 正向 Top 8
kegg.up <- kegg.sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 8)

# 负向 Top 8
kegg.down <- kegg.sig %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 8)

kegg.plot <- bind_rows(
  kegg.up,
  kegg.down
) %>%
  arrange(NES)

ggplot(
  kegg.plot,
  aes(
    x = NES,
    y = reorder(pathway, NES)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    aes(
      size = -log10(p.adjust),
      color = NES
    ),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      margin = margin(t = 8)
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = "NES",
    size = "-log10 adjusted P"
  )







library(dplyr)

rank.day <- partial.cor.results.iav271812 %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(rho_day_adjusted),
    is.finite(rho_day_adjusted)
  ) %>%
  arrange(desc(rho_day_adjusted))

gsea.rank.day <- rank.day$rho_day_adjusted
names(gsea.rank.day) <- rank.day$gene

gsea.rank.day <- gsea.rank.day[
  !duplicated(names(gsea.rank.day))
]

gsea.rank.day <- sort(
  gsea.rank.day,
  decreasing = TRUE
)



head(gsea.rank.day, 20)
tail(gsea.rank.day, 20)
length(gsea.rank.day)





rank.raw <- cor.results.iav271812 %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(rho),
    is.finite(rho)
  ) %>%
  arrange(desc(rho))

gsea.rank.raw <- rank.raw$rho

names(gsea.rank.raw) <- rank.raw$gene

gsea.rank.raw <- gsea.rank.raw[
  !duplicated(names(gsea.rank.raw))
]

gsea.rank.raw <- sort(
  gsea.rank.raw,
  decreasing = TRUE
)




library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

kegg.map <- bitr(
  names(gsea.rank.raw),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)



head(kegg.map)
nrow(kegg.map)


kegg.rank.raw.df <- data.frame(
  SYMBOL = names(gsea.rank.raw),
  score = as.numeric(gsea.rank.raw),
  stringsAsFactors = FALSE
)

kegg.rank.raw.df <- kegg.rank.raw.df %>%
  inner_join(
    kegg.map,
    by = "SYMBOL"
  )


kegg.rank.raw <- kegg.rank.raw.df$score

names(kegg.rank.raw) <-
  kegg.rank.raw.df$ENTREZID

kegg.rank.raw <- sort(
  kegg.rank.raw,
  decreasing = TRUE
)



gsea.kegg.raw <- gseKEGG(
  geneList = kegg.rank.raw,
  organism = "mmu",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)





kegg.raw.df <- as.data.frame(
  gsea.kegg.raw
)



kegg.raw.df %>%
  arrange(p.adjust) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  ) %>%
  head(30)



kegg.raw.df %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  ) %>%
  head(30)






kegg.raw.sig <- kegg.raw.df %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  dplyr::arrange(p.adjust)

nrow(kegg.raw.sig)


kegg.raw.sig %>%
  dplyr::select(
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  )



kegg.raw.sig %>%
  dplyr::mutate(
    direction = ifelse(
      NES > 0,
      "Positive",
      "Negative"
    )
  ) %>%
  dplyr::count(direction)


kegg.raw.sig %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(dplyr::desc(NES)) %>%
  dplyr::select(
    Description,
    NES,
    p.adjust,
    setSize
  )




kegg.raw.sig %>%
  dplyr::filter(NES < 0) %>%
  dplyr::arrange(NES) %>%
  dplyr::select(
    Description,
    NES,
    p.adjust,
    setSize
  )



library(ggplot2)
library(stringr)

kegg.up <- kegg.raw.sig %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(dplyr::desc(NES)) %>%
  dplyr::slice_head(n = 10)

kegg.down <- kegg.raw.sig %>%
  dplyr::filter(NES < 0) %>%
  dplyr::arrange(NES) %>%
  dplyr::slice_head(n = 8)

kegg.plot <- dplyr::bind_rows(
  kegg.down,
  kegg.up
) %>%
  dplyr::mutate(
    pathway = stringr::str_wrap(
      Description,
      width = 38
    )
  )





ggplot(
  kegg.plot,
  aes(
    x = NES,
    y = reorder(pathway, NES)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    aes(
      size = -log10(p.adjust),
      color = NES
    ),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "grey90",
    high = "#D73027",
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = "NES",
    size = "-log10 adjusted P"
  )


kegg.raw.df %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  ) %>%
  head(30)



kegg.raw.df %>%
  dplyr::filter(
    grepl(
      "JAK|STAT|Toll|NF-kappa|NF-kB|PI3K|Akt|MAPK|complement|coagulation|interferon|immune|inflamm|cytokine",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  ) %>%
  dplyr::arrange(p.adjust)





kegg.raw.df %>%
  dplyr::filter(
    grepl(
      "MAPK",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust
  )



library(enrichplot)
library(ggplot2)

kegg.ids <- kegg.raw.df %>%
  dplyr::filter(
    Description %in% c(
      "NF-kappa B signaling pathway",
      "Toll-like receptor signaling pathway",
      "JAK-STAT signaling pathway",
      "PI3K-Akt signaling pathway"
    )
  ) %>%
  dplyr::pull(ID)

kegg.ids


p.kegg <- gseaplot2(
  gsea.kegg.raw,
  geneSetID = kegg.ids,
  title = NULL,
  color = c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488"
  ),
  base_size = 11,
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)

p.kegg



kegg.raw.df %>%
  dplyr::filter(ID %in% kegg.ids) %>%
  dplyr::select(
    ID,
    Description,
    NES,
    p.adjust
  )


head(gsea.kegg.raw@geneList)
length(gsea.kegg.raw@geneList)
head(names(gsea.kegg.raw@geneList))



partial.cor.results.iav271812


# Historical interactive inspection of .Rhistory omitted from upload copy.
# This is not a step required for the published analyses.



day.factor <- factor(meta.iav271812$day)

ism1.resid.factor <- residuals(
  lm(
    ism1.expression.iav271812 ~ day.factor
  )
)

ehd3.expression <- as.numeric(
  vst.matrix.iav271812["Ehd3", ]
)

ehd3.resid.factor <- residuals(
  lm(
    ehd3.expression ~ day.factor
  )
)

cor.test(
  ism1.resid.factor,
  ehd3.resid.factor,
  method = "spearman",
  exact = FALSE
)






############################################################
## B1. Day-adjusted ISM1-associated correlations
## Day treated as categorical variable
############################################################

day.factor.iav271812 <- factor(
  meta.iav271812$day
)

# 检查 Day
table(day.factor.iav271812)


ism1.residual.factorDay.iav271812 <- residuals(
  lm(
    ism1.expression.iav271812 ~
      day.factor.iav271812
  )
)



############################################################
## B3. Calculate Day-adjusted Spearman correlations
############################################################

genes.iav271812 <- setdiff(
  rownames(vst.matrix.iav271812),
  "Ism1"
)

partial.cor.factorDay.iav271812 <- lapply(
  
  genes.iav271812,
  
  function(gene) {
    
    gene.expression <- as.numeric(
      vst.matrix.iav271812[gene, ]
    )
    
    # Remove Day effect
    gene.residual <- residuals(
      lm(
        gene.expression ~
          day.factor.iav271812
      )
    )
    
    # Spearman correlation between
    # Day-adjusted ISM1 and Day-adjusted gene
    test <- suppressWarnings(
      cor.test(
        ism1.residual.factorDay.iav271812,
        gene.residual,
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
      gene = gene,
      rho_day_adjusted = unname(test$estimate),
      pvalue_day_adjusted = test$p.value,
      stringsAsFactors = FALSE
    )
  }
)

partial.cor.factorDay.iav271812 <-
  dplyr::bind_rows(
    partial.cor.factorDay.iav271812
  )




partial.cor.factorDay.iav271812$padj_day_adjusted <-
  p.adjust(
    partial.cor.factorDay.iav271812$pvalue_day_adjusted,
    method = "BH"
  )



partial.cor.factorDay.iav271812 <-
  partial.cor.factorDay.iav271812 |>
  dplyr::arrange(
    desc(rho_day_adjusted)
  )


partial.cor.factorDay.iav271812 |>
  dplyr::filter(
    rho_day_adjusted > 0
  ) |>
  dplyr::arrange(
    padj_day_adjusted
  ) |>
  head(30)

partial.cor.factorDay.iav271812 |>
  dplyr::filter(
    rho_day_adjusted < 0
  ) |>
  dplyr::arrange(
    padj_day_adjusted
  ) |>
  head(30)


partial.cor.factorDay.iav271812 |>
  dplyr::filter(
    gene == "Ehd3"
  )





length(geneList.day.adjusted.iav271812)

head(geneList.day.adjusted.iav271812)

tail(geneList.day.adjusted.iav271812)

sum(is.na(geneList.day.adjusted.iav271812))

sum(!is.finite(geneList.day.adjusted.iav271812))



write.csv(
  partial.cor.factorDay.iav271812,
  file = file.path(
    result.dir.271812,
    "GSE271812_Ism1_IAV_Spearman_factorDay_adjusted.csv"
  ),
  row.names = FALSE
)





geneList.day.adjusted.iav271812 <-
  partial.cor.factorDay.iav271812$rho_day_adjusted

names(geneList.day.adjusted.iav271812) <-
  partial.cor.factorDay.iav271812$gene


head(geneList.day.adjusted.iav271812)


sum(duplicated(names(geneList.day.adjusted.iav271812)))


geneList.day.adjusted.iav271812 <-
  sort(
    geneList.day.adjusted.iav271812,
    decreasing = TRUE
  )


head(
  geneList.day.adjusted.iav271812,
  20
)

length(
  geneList.day.adjusted.iav271812
)
is.unsorted(
  geneList.day.adjusted.iav271812,
  strictly = FALSE
)



nrow(partial.cor.factorDay.iav271812)

length(rownames(vst.matrix.iav271812))

length(setdiff(
  rownames(vst.matrix.iav271812),
  "Ism1"
))


geneList.day.adjusted.iav271812 <-
  sort(
    geneList.day.adjusted.iav271812,
    decreasing = TRUE
  )
is.unsorted(
  geneList.day.adjusted.iav271812,
  strictly = FALSE
)
tail(
  geneList.day.adjusted.iav271812,
  20
)
head(
  geneList.day.adjusted.iav271812,
  20
)
all(
  diff(geneList.day.adjusted.iav271812) <= 0
)
sum(is.na(geneList.day.adjusted.iav271812))
sum(!is.finite(geneList.day.adjusted.iav271812))



length(gsea.kegg.raw@geneList)

head(names(gsea.kegg.raw@geneList))

tail(names(gsea.kegg.raw@geneList))

head(names(geneList.day.adjusted.iav271812))

sum(
  names(geneList.day.adjusted.iav271812) %in%
    names(gsea.kegg.raw@geneList)
)

sum(
  names(gsea.kegg.raw@geneList) %in%
    names(geneList.day.adjusted.iav271812)
)



class(gsea.kegg.raw)

gsea.kegg.raw@params

head(kegg.geneSets)

str(kegg.geneSets)






slotNames(gsea.kegg.raw)









############################################################
## B16. Build Day-adjusted Entrez ranking
############################################################

day.rank.df.iav271812 <-
  data.frame(
    SYMBOL = names(geneList.day.adjusted.iav271812),
    rho_day_adjusted =
      as.numeric(geneList.day.adjusted.iav271812),
    stringsAsFactors = FALSE
  ) |>
  dplyr::left_join(
    id.map.day.iav271812,
    by = "SYMBOL"
  )


head(day.rank.df.iav271812)

head(
  gsea.kegg.raw@result
)

head(
  gsea.kegg.raw@geneSets
)



############################################################
## B16. Build Day-adjusted Entrez ranking
############################################################

day.rank.df.iav271812 <-
  data.frame(
    SYMBOL = names(geneList.day.adjusted.iav271812),
    rho_day_adjusted =
      as.numeric(geneList.day.adjusted.iav271812),
    stringsAsFactors = FALSE
  ) |>
  dplyr::left_join(
    id.map.day.iav271812,
    by = "SYMBOL"
  )



head(day.rank.df.iav271812)



day.rank.df.iav271812 <-
  day.rank.df.iav271812 |>
  dplyr::filter(
    !is.na(ENTREZID)
  )


nrow(day.rank.df.iav271812)

sum(
  duplicated(
    day.rank.df.iav271812$ENTREZID
  )
)


geneList.day.adjusted.entrez.iav271812 <-
  day.rank.df.iav271812$rho_day_adjusted

names(geneList.day.adjusted.entrez.iav271812) <-
  day.rank.df.iav271812$ENTREZID



geneList.day.adjusted.entrez.iav271812 <-
  sort(
    geneList.day.adjusted.entrez.iav271812,
    decreasing = TRUE
  )


length(
  geneList.day.adjusted.entrez.iav271812
)



sum(
  names(geneList.day.adjusted.entrez.iav271812) %in%
    names(gsea.kegg.raw@geneList)
)
sum(
  names(gsea.kegg.raw@geneList) %in%
    names(geneList.day.adjusted.entrez.iav271812)
)




geneList.day.adjusted.entrez.iav271812["57440"]

id.map.day.iav271812 |>
  dplyr::filter(
    SYMBOL == "Ehd3"
  )




length(geneList.day.adjusted.entrez.iav271812)

sum(
  names(geneList.day.adjusted.entrez.iav271812) %in%
    names(gsea.kegg.raw@geneList)
)

sum(
  names(gsea.kegg.raw@geneList) %in%
    names(geneList.day.adjusted.entrez.iav271812)
)

all(
  diff(geneList.day.adjusted.entrez.iav271812) <= 0
)



kegg.term2gene.iav271812 <-
  do.call(
    rbind,
    lapply(
      names(gsea.kegg.raw@geneSets),
      function(pathway) {
        
        data.frame(
          term = pathway,
          gene = gsea.kegg.raw@geneSets[[pathway]],
          stringsAsFactors = FALSE
        )
        
      }
    )
  )



head(kegg.term2gene.iav271812)
dim(kegg.term2gene.iav271812)


library(clusterProfiler)

gsea.kegg.day.adjusted.iav271812 <-
  GSEA(
    geneList =
      geneList.day.adjusted.entrez.iav271812,
    
    TERM2GENE =
      kegg.term2gene.iav271812,
    
    pvalueCutoff = 0.05,
    
    pAdjustMethod = "BH",
    
    exponent = 1,
    
    minGSSize = 10,
    
    maxGSSize = 500,
    
    eps = 1e-10
  )




kegg.day.adjusted.df <-
  as.data.frame(
    gsea.kegg.day.adjusted.iav271812
  )


kegg.day.adjusted.df |>
  dplyr::arrange(p.adjust) |>
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  ) |>
  head(30)





kegg.day.adjusted.df |>
  dplyr::filter(
    ID %in% c(
      "mmu04064",
      "mmu04620",
      "mmu04630",
      "mmu04151"
    )
  ) |>
  dplyr::select(
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    setSize
  )



raw.kegg.compare <-
  kegg.raw.df |>
  dplyr::filter(
    ID %in% c(
      "mmu04064",
      "mmu04620",
      "mmu04630",
      "mmu04151"
    )
  ) |>
  dplyr::select(
    ID,
    Description,
    NES.raw = NES,
    FDR.raw = p.adjust
  )

day.kegg.compare <-
  kegg.day.adjusted.df |>
  dplyr::filter(
    ID %in% c(
      "mmu04064",
      "mmu04620",
      "mmu04630",
      "mmu04151"
    )
  ) |>
  dplyr::select(
    ID,
    Description,
    NES.day = NES,
    FDR.day = p.adjust
  )

kegg.compare.iav271812 <-
  dplyr::full_join(
    raw.kegg.compare,
    day.kegg.compare,
    by = c("ID", "Description")
  )

kegg.compare.iav271812






raw.kegg.compare <-
  kegg.raw.df |>
  dplyr::filter(
    ID %in% c(
      "mmu04064",
      "mmu04620",
      "mmu04630",
      "mmu04151"
    )
  ) |>
  dplyr::select(
    ID,
    Description,
    NES.raw = NES,
    FDR.raw = p.adjust
  )

day.kegg.compare <-
  kegg.day.adjusted.df |>
  dplyr::filter(
    ID %in% c(
      "mmu04064",
      "mmu04620",
      "mmu04630",
      "mmu04151"
    )
  ) |>
  dplyr::select(
    ID,
    NES.day = NES,
    FDR.day = p.adjust
  )

kegg.compare.iav271812 <-
  dplyr::full_join(
    raw.kegg.compare,
    day.kegg.compare,
    by = "ID"
  )

kegg.compare.iav271812




############################################################
## Figure 3A
## Raw ISM1-associated GO-BP GSEA
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

go.rank.raw <- cor.results.iav271812 %>%
  dplyr::filter(
    !is.na(gene),
    gene != "",
    !is.na(rho),
    is.finite(rho)
  )

geneList.go.raw <- go.rank.raw$rho
names(geneList.go.raw) <- go.rank.raw$gene

geneList.go.raw <- geneList.go.raw[
  !duplicated(names(geneList.go.raw))
]

geneList.go.raw <- sort(
  geneList.go.raw,
  decreasing = TRUE
)





gsea.go.raw.iav271812 <- gseGO(
  geneList = geneList.go.raw,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 1e-10,
  verbose = FALSE
)











library(ggplot2)
library(dplyr)
library(patchwork)




pathways.figure3DG <- c(
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION"
)

stopifnot(
  all(
    pathways.figure3DG %in%
      colnames(gsva.plot.data.271812)
  )
)







make.gsva.correlation.plot <- function(
    data,
    pathway,
    title,
    panel.label
) {
  
  x <- data$Ism1_VST
  y <- data[[pathway]]
  
  test <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(test$estimate)
  pval <- test$p.value
  
  stat.label <- paste0(
    "\u03c1 = ",
    sprintf("%.2f", rho),
    ", P = ",
    ifelse(
      pval < 0.001,
      format.pval(pval, digits = 2, eps = 0.001),
      sprintf("%.3f", pval)
    )
  )
  
  ggplot(
    data,
    aes(
      x = Ism1_VST,
      y = .data[[pathway]]
    )
  ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8,
      color = "#365A8C",
      fill = "#C9D3E6",
      alpha = 0.55
    ) +
    
    geom_point(
      size = 3.1,
      shape = 21,
      stroke = 0.4,
      fill = "#E97868",
      color = "white"
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat.label,
      hjust = -0.05,
      vjust = 1.4,
      size = 4.4,
      fontface = "bold"
    ) +
    
    labs(
      title = title,
      x = "ISM1 expression (VST)",
      y = "GSVA score",
      tag = panel.label
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      
      plot.tag = element_text(
        face = "bold",
        size = 18
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 11.5
      ),
      
      axis.text = element_text(
        color = "black",
        size = 10
      ),
      
      axis.line = element_line(
        linewidth = 0.7
      ),
      
      plot.margin = margin(
        10, 12, 10, 10
      )
    )
}







p.D <- make.gsva.correlation.plot(
  data = gsva.plot.data.271812,
  pathway =
    "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  title =
    "TOLL-LIKE RECEPTOR\nSIGNALING PATHWAY",
  panel.label = "D"
)








library(ggplot2)
library(dplyr)
library(patchwork)

make.gsva.correlation.plot.pretty <- function(
    data,
    pathway,
    title,
    panel.label
) {
  
  x <- data$Ism1_VST
  y <- data[[pathway]]
  
  test <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(test$estimate)
  pval <- test$p.value
  
  stat.label <- paste0(
    "Spearman \u03c1 = ",
    sprintf("%.2f", rho),
    "\nP = ",
    ifelse(
      pval < 0.001,
      "< 0.001",
      sprintf("%.3f", pval)
    )
  )
  
  ggplot(
    data,
    aes(
      x = Ism1_VST,
      y = .data[[pathway]]
    )
  ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "#2F5D8C",
      fill = "#DCE6F2",
      alpha = 0.7
    ) +
    
    geom_point(
      size = 3.4,
      shape = 21,
      stroke = 0.5,
      fill = "#E76F51",
      color = "white"
    ) +
    
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = stat.label,
      hjust = -0.05,
      vjust = 1.15,
      size = 3.9,
      fontface = "bold",
      label.size = 0,
      fill = scales::alpha("white", 0.8)
    ) +
    
    labs(
      title = title,
      x = "ISM1 expression (VST)",
      y = "GSVA score",
      tag = panel.label
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      
      plot.tag = element_text(
        face = "bold",
        size = 17
      ),
      
      axis.title.x = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(t = 8)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(r = 8)
      ),
      
      axis.text = element_text(
        color = "black",
        size = 9.5
      ),
      
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      
      plot.margin = margin(
        10, 12, 10, 12
      )
    )
}




p.D <- make.gsva.correlation.plot.pretty(
  data = gsva.plot.data.271812,
  pathway = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  title = "Toll-like receptor\nsignaling pathway",
  panel.label = "D"
)

p.E <- make.gsva.correlation.plot.pretty(
  data = gsva.plot.data.271812,
  pathway = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  title = "Complement and\ncoagulation cascades",
  panel.label = "E"
)

p.F <- make.gsva.correlation.plot.pretty(
  data = gsva.plot.data.271812,
  pathway = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  title = "JAK-STAT\nsignaling pathway",
  panel.label = "F"
)

p.G <- make.gsva.correlation.plot.pretty(
  data = gsva.plot.data.271812,
  pathway = "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
  title = "Positive regulation of type I\ninterferon production",
  panel.label = "G"
)








p.figure3DG.pretty <-
  (p.D | p.E) /
  (p.F | p.G) +
  
  plot_annotation(
    title = "ISM1 expression is associated with antiviral and inflammatory pathway activity"
  ) &
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

p.figure3DG.pretty




ggsave(
  "GSE271812_ISM1_GSVA_Figure3DG_pretty.png",
  p.figure3DG.pretty,
  width = 9,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  "GSE271812_ISM1_GSVA_Figure3DG_pretty.pdf",
  p.figure3DG.pretty,
  width = 9,
  height = 8,
  bg = "white"
)







library(dplyr)
library(ggplot2)
library(patchwork)

############################################################
## Day 4 only: Mock + IAV
############################################################

meta.day4.271812 <- meta.gsva.271812 |>
  dplyr::filter(
    day == 4,
    infection_group %in% c("Mock", "IAV")
  )

meta.day4.271812 <- as.data.frame(meta.day4.271812)

print(
  meta.day4.271812[, c(
    "sample",
    "day",
    "infection_group"
  )]
)

table(meta.day4.271812$infection_group)









coldata.day4.gsva <- data.frame(
  condition = metadata.primary271812$condition,
  row.names = metadata.primary271812$sample
)

count.day4.gsva <- count.matrix.primary271812[
  ,
  metadata.primary271812$sample,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(count.day4.gsva),
    rownames(coldata.day4.gsva)
  )
)

dds.day4.gsva <- DESeq2::DESeqDataSetFromMatrix(
  countData = count.day4.gsva,
  colData = coldata.day4.gsva,
  design = ~ condition
)






keep.day4.gsva <- rowSums(
  DESeq2::counts(dds.day4.gsva) >= 5
) >= 3

dds.day4.gsva <- dds.day4.gsva[
  keep.day4.gsva,
]

cat(
  "Genes retained:",
  nrow(dds.day4.gsva),
  "\n"
)






vsd.day4.gsva <- DESeq2::vst(
  dds.day4.gsva,
  blind = TRUE
)

vst.matrix.day4.gsva <- SummarizedExperiment::assay(
  vsd.day4.gsva
)

dim(vst.matrix.day4.gsva)






colnames(vst.matrix.day4.gsva)

metadata.primary271812[, c(
  "sample",
  "condition",
  "replicate"
)]


library(msigdbr)
library(dplyr)

msig.kegg.day4 <- msigdbr::msigdbr(
  species = "Mus musculus",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
)

msig.go.day4 <- msigdbr::msigdbr(
  species = "Mus musculus",
  collection = "C5",
  subcollection = "GO:BP"
)






target.pathways.day4 <- c(
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION"
)




target.genes.day4 <- dplyr::bind_rows(
  
  msig.kegg.day4 |>
    dplyr::filter(
      gs_name %in% target.pathways.day4
    ) |>
    dplyr::select(
      gs_name,
      gene_symbol
    ),
  
  msig.go.day4 |>
    dplyr::filter(
      gs_name %in% target.pathways.day4
    ) |>
    dplyr::select(
      gs_name,
      gene_symbol
    )
  
) |>
  dplyr::distinct()





gene.set.list.day4 <- split(
  target.genes.day4$gene_symbol,
  target.genes.day4$gs_name
)




gene.set.list.day4 <- lapply(
  gene.set.list.day4,
  function(x) {
    intersect(
      x,
      rownames(vst.matrix.day4.gsva)
    )
  }
)

sapply(
  gene.set.list.day4,
  length
)




gene.set.list.day4 <- gene.set.list.day4[
  sapply(
    gene.set.list.day4,
    length
  ) >= 10
]




library(GSVA)

gsva.param.day4 <- GSVA::gsvaParam(
  exprData = as.matrix(
    vst.matrix.day4.gsva
  ),
  geneSets = gene.set.list.day4,
  kcdf = "Gaussian"
)

gsva.result.day4 <- GSVA::gsva(
  gsva.param.day4,
  verbose = TRUE
)

gsva.matrix.day4 <- as.matrix(
  gsva.result.day4
)




dim(gsva.matrix.day4)

rownames(gsva.matrix.day4)

colnames(gsva.matrix.day4)



gsva.day4.plot.data <- data.frame(
  
  sample =
    colnames(
      gsva.matrix.day4
    ),
  
  ISM1 =
    as.numeric(
      vst.matrix.day4.gsva[
        "Ism1",
        colnames(
          gsva.matrix.day4
        )
      ]
    ),
  
  stringsAsFactors = FALSE
  
) |>
  dplyr::left_join(
    
    metadata.primary271812 |>
      dplyr::select(
        sample,
        condition,
        replicate
      ),
    
    by = "sample"
  )






for (
  p in rownames(gsva.matrix.day4)
) {
  
  gsva.day4.plot.data[[p]] <-
    as.numeric(
      gsva.matrix.day4[
        p,
        gsva.day4.plot.data$sample
      ]
    )
}





table(
  gsva.day4.plot.data$condition
)




make.day4.gsva.plot <- function(
    data,
    pathway,
    title,
    panel
) {
  
  test <- cor.test(
    data$ISM1,
    data[[pathway]],
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(
    test$estimate
  )
  
  pval <- test$p.value
  
  stat.text <- paste0(
    "\u03c1 = ",
    sprintf("%.2f", rho),
    ", P = ",
    ifelse(
      pval < 0.001,
      "<0.001",
      sprintf("%.3f", pval)
    )
  )
  
  ggplot(
    data,
    aes(
      x = ISM1,
      y = .data[[pathway]]
    )
  ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "#355C8A",
      fill = "#D7DFEC",
      alpha = 0.7
    ) +
    
    geom_point(
      aes(
        fill = condition
      ),
      shape = 21,
      size = 3.5,
      stroke = 0.5,
      color = "white"
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "#B8B8B8",
        "Influenza" = "#E87563"
      )
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat.text,
      hjust = -0.03,
      vjust = 1.35,
      fontface = "bold",
      size = 4
    ) +
    
    labs(
      title = title,
      tag = panel,
      x = "ISM1 expression (VST)",
      y = "GSVA score",
      fill = NULL
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 11
      ),
      plot.tag = element_text(
        face = "bold",
        size = 17
      ),
      axis.title = element_text(
        face = "bold"
      ),
      legend.position = "top"
    )
}



table(gsva.day4.plot.data$condition)



library(ggplot2)
library(dplyr)
library(patchwork)

make.day4.gsva.plot.pretty <- function(
    data,
    pathway,
    title,
    panel
) {
  
  # 检查 pathway 是否存在
  if (!pathway %in% colnames(data)) {
    stop(
      paste0(
        "Pathway not found in data: ",
        pathway
      )
    )
  }
  
  # Spearman correlation
  test <- cor.test(
    data$ISM1,
    data[[pathway]],
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(test$estimate)
  pval <- test$p.value
  
  stat.text <- paste0(
    "\u03c1 = ",
    sprintf("%.2f", rho),
    ", P = ",
    ifelse(
      pval < 0.001,
      "<0.001",
      sprintf("%.3f", pval)
    )
  )
  
  # plot
  ggplot(
    data,
    aes(
      x = ISM1,
      y = .data[[pathway]]
    )
  ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "#355C8A",
      fill = "#D7DFEC",
      alpha = 0.70
    ) +
    
    geom_point(
      aes(
        fill = condition
      ),
      shape = 21,
      size = 3.8,
      stroke = 0.6,
      color = "white"
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "#9E9E9E",
        "Influenza" = "#E76F51"
      )
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat.text,
      hjust = -0.03,
      vjust = 1.35,
      size = 4.2,
      fontface = "bold"
    ) +
    
    labs(
      title = title,
      tag = panel,
      x = "ISM1 expression (VST)",
      y = "GSVA score",
      fill = NULL
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        lineheight = 0.95,
        margin = margin(b = 7)
      ),
      
      plot.tag = element_text(
        face = "bold",
        size = 18
      ),
      
      axis.title.x = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(t = 8)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(r = 8)
      ),
      
      axis.text = element_text(
        size = 9.5,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      
      legend.position = "top",
      
      legend.text = element_text(
        size = 10
      ),
      
      plot.margin = margin(
        12, 14, 12, 14
      )
    )
}




p.D.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  title = "Toll-like receptor\nsignaling pathway",
  panel = "D"
)

p.E.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  title = "Complement and\ncoagulation cascades",
  panel = "E"
)

p.F.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  title = "JAK-STAT\nsignaling pathway",
  panel = "F"
)

p.G.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway =
    "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
  title =
    "Positive regulation of type I\ninterferon production",
  panel = "G"
)






p.figure3DG.day4.pretty <-
  (p.D.day4 | p.E.day4) /
  (p.F.day4 | p.G.day4) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

p.figure3DG.day4.pretty





ggsave(
  filename = file.path(
    figure.dir.271812,
    "GSE271812_ISM1_D4_GSVA_correlation_4pathways.png"
  ),
  plot = p.figure3DG.day4.pretty,
  width = 9,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure.dir.271812,
    "GSE271812_ISM1_D4_GSVA_correlation_4pathways.pdf"
  ),
  plot = p.figure3DG.day4.pretty,
  width = 9,
  height = 8,
  bg = "white"
)










library(ggplot2)
library(patchwork)

make.day4.gsva.plot.pretty <- function(
    data,
    pathway,
    title
) {
  
  test <- cor.test(
    data$ISM1,
    data[[pathway]],
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(test$estimate)
  pval <- test$p.value
  
  stat.text <- paste0(
    "\u03c1 = ",
    sprintf("%.2f", rho),
    ", P = ",
    ifelse(
      pval < 0.001,
      "<0.001",
      sprintf("%.3f", pval)
    )
  )
  
  ggplot(
    data,
    aes(
      x = ISM1,
      y = .data[[pathway]]
    )
  ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "#355C8A",
      fill = "#D7DFEC",
      alpha = 0.70
    ) +
    
    geom_point(
      aes(fill = condition),
      shape = 21,
      size = 4,
      stroke = 0.6,
      color = "white"
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "#A6A6A6",
        "Influenza" = "#E97766"
      )
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat.text,
      hjust = -0.03,
      vjust = 1.35,
      size = 4.2,
      fontface = "bold"
    ) +
    
    labs(
      title = title,
      x = "ISM1 expression (VST)",
      y = "GSVA score",
      fill = NULL
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        lineheight = 0.95,
        margin = margin(b = 8)
      ),
      
      axis.title.x = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(t = 8)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(r = 8)
      ),
      
      axis.text = element_text(
        size = 9.5,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      
      legend.position = "top",
      
      legend.text = element_text(
        size = 10
      ),
      
      plot.margin = margin(
        12, 14, 12, 14
      )
    )
}







p.TLR.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  title = "Toll-like receptor\nsignaling pathway"
)

p.Coag.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  title = "Complement and\ncoagulation cascades"
)

p.JAK.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  title = "JAK-STAT\nsignaling pathway"
)

p.IFN.day4 <- make.day4.gsva.plot.pretty(
  data = gsva.day4.plot.data,
  pathway =
    "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
  title =
    "Positive regulation of type I\ninterferon production"
)







p.figure.day4 <-
  (p.TLR.day4 | p.Coag.day4) /
  (p.JAK.day4 | p.IFN.day4) +
  patchwork::plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

p.figure.day4





ggsave(
  file.path(
    figure.dir.271812,
    "GSE271812_ISM1_D4_GSVA_correlation_4pathways.png"
  ),
  p.figure.day4,
  width = 9,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure.dir.271812,
    "GSE271812_ISM1_D4_GSVA_correlation_4pathways.pdf"
  ),
  p.figure.day4,
  width = 9,
  height = 8,
  bg = "white"
)









############################################################
## Remove D / E / F / G labels
############################################################

p.D.clean <- p.D +
  labs(tag = NULL) +
  theme(plot.tag = element_blank())

p.E.clean <- p.E +
  labs(tag = NULL) +
  theme(plot.tag = element_blank())

p.F.clean <- p.F +
  labs(tag = NULL) +
  theme(plot.tag = element_blank())

p.G.clean <- p.G +
  labs(tag = NULL) +
  theme(plot.tag = element_blank())




library(patchwork)

p.figure3DG.clean <-
  (p.D.clean | p.E.clean) /
  (p.F.clean | p.G.clean) +
  
  plot_annotation(
    title =
      "ISM1 expression is associated with antiviral and inflammatory pathway activity"
  ) &
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    )
  )

p.figure3DG.clean




ggsave(
  filename = file.path(
    figure.dir.271812,
    "GSE271812_ISM1_GSVA_correlation_4pathways_clean.png"
  ),
  plot = p.figure3DG.clean,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure.dir.271812,
    "GSE271812_ISM1_GSVA_correlation_4pathways_clean.pdf"
  ),
  plot = p.figure3DG.clean,
  width = 10,
  height = 8,
  bg = "white"
)