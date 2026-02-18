args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript diff_expression_swe.R <raw_count_matrix_file> <col_data_file>")
}
matrix_file <- args[1]
col_data_file <- args[2]

library("DESeq2")
library("ggplot2")
library(ggrepel)
library(EnhancedVolcano)
library("pheatmap")

# Load data
counts_data = read.table(matrix_file, header=T, row.names=1, com='');
print("Loading raw count matrix: Done")
col_data <- read.table(
  col_data_file,
  header = TRUE,
  row.names = 1,
  com = "",
  stringsAsFactors = FALSE,
  na.strings = ""
)
col_data$temperature <- factor(col_data$temperature)
col_data$population  <- factor(col_data$population)
col_data$family      <- factor(col_data$family)
col_data$baltic_side <- factor(col_data$baltic_side)
col_data$region <- factor(col_data$region)
print("Loading col data: Done")

# Reorder col_data rows to match counts_data column order
col_data <- col_data[colnames(counts_data), ]
col_all <- col_data

# Subset count matrix to the same samples
cts <- counts_data[, rownames(col_all)]
# Round matrix
cts <- round(cts)
# Check matrix and coldata structure
cat("Check matrix and coldata structure. Is it correct?", ncol(cts) == nrow(col_all), "\n")

# Construct DESeq2 dataset
# Put the variable of interest at the end of the formula and make sure the control level is the first level.
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = col_all,
  design = ~ temperature + family + region
)

# Pre filtering
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# Set Factor Level
dds$region <- relevel(dds$region, ref = "South")

# Run DESeq
print("Running DESeq")
dds <- DESeq(dds)
dds_EN <- dds
dds_EN$region <- relevel(dds_EN$region, ref = "North")
dds_EN <- DESeq(dds_EN)
print("Obtaining results")
res_East_South <- results(dds, name = "region_East_vs_South")
res_North_South <- results(dds, name = "region_North_vs_South")
res_East_North <- results(dds_EN, name = "region_East_vs_North")
print("Running DESeq: DONE")

print("Calculating Shrinkage of effect size")
resLFC_East_South  <- lfcShrink(dds, coef = "region_East_vs_South",  type = "apeglm")
resLFC_North_South <- lfcShrink(dds, coef = "region_North_vs_South", type = "apeglm")
resLFC_East_North <- lfcShrink(dds_EN, coef = "region_East_vs_North", type = "apeglm")

print("Saving Results")
save(list=ls(all=TRUE), file="05_Region.RData")
print("Saving Results: DONE")

print("-----------------------------------------------------------------")
warnings()

summary(res_East_South)
summary(res_North_South)
summary(res_East_North)
summary(resLFC_East_South)
summary(resLFC_North_South)
summary(resLFC_East_North)


# Convert to data frame
res_df_ES <- as.data.frame(resLFC_East_South)
res_df_NS <- as.data.frame(resLFC_North_South)
res_df_EN <- as.data.frame(resLFC_East_North)

# ------------------------------------------------------------------------------
# MA PLOT

# Define categories
res_df$diffexp_sig <- "Not significant"

sig_idx <- !is.na(res_df$padj) & res_df$padj < 0.05

res_df$diffexp_sig[sig_idx & res_df$log2FoldChange >= 1]  <- "Upregulated"
res_df$diffexp_sig[sig_idx & res_df$log2FoldChange <= -1] <- "Downregulated"
res_df$diffexp_sig[sig_idx & abs(res_df$log2FoldChange) < 1] <- "Significant, small effect"

# MA plot
res_df$diffexp_sig <- factor(
  res_df$diffexp_sig,
  levels = c(
    "Downregulated",
    "Upregulated",
    "Significant, small effect",
    "Not significant"
  )
)

res_df_plot <- res_df[!is.na(res_df$padj), ]

ggplot(res_df_plot, aes(x = baseMean, y = log2FoldChange, color = diffexp_sig)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Downregulated" = "blue",
      "Upregulated" = "red",
      "Significant, small effect" = "gray40",
      "Not significant" = "gray80"
    ),
    name = "Differential expression"
  ) +
  geom_hline(yintercept = 0) +
  labs(
    title = "MA plot: Region effect (East vs North)",
    x = "Mean of normalized counts",
    y = "Log2 fold change (shrunken)"
  ) +
  coord_cartesian(ylim = c(-20, 20)) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Volcano Plot

EnhancedVolcano(
  resLFC_East_North,
  lab = NA,  # no labels by default
  x = 'log2FoldChange',
  y = 'padj',

  xlim = c(-16,16),
  ylim = c(0, 20),
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  pointSize = 1.8,
  labSize = 3.0,
  
  title = 'Volcano plot: region effect (East vs North)',
  subtitle = 'padj < 0.05, |log2FC| ≥ 1',
  
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4,
  
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)

# ------------------------------------------------------------------------------
# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
# ------------------------------------------------------------------------------

# PCA

pop_colors <- c(
  # Kalmar (South)
  "Ka_15"  = "#FFDE52",
  "Ka_20"  = "#E8B802",
  # Uppsala (South)
  "Upp_15" = "#FF7070",
  "Upp_20" = "#BF1515",
  # Latvia (South)
  "L_15"  = "#00CDFF",
  "L_20"  = "#00718A",
  # Estonia (South)
  "E_15" = "#B7FF5E",
  "E_20" = "#488200",
  # Finland (North)
  "C_Fin_15" = "#FF33E4",  
  "C_Fin_20" = "#7D006A",
  # North populations
  "NA_15"  = "#78B1FF",
  "NA_20"  = "#0052BF",
  "NL_15"  = "#59FFC0",
  "NL_20"  = "#00B06C",
  "VF_15"  = "#bcbddc",
  "VF_20"  = "#54278f"
)

get_pca_data <- function(vsd, pcs) {
  pcaData <- plotPCA(
    vsd,
    intgroup = c("population", "temperature", "region"),
    returnData = TRUE,
    pcsToUse = pcs
  )
  
  pcaData$pop_temp <- factor(
    paste(pcaData$population, pcaData$temperature, sep = "_"),
    levels = c(
      "L_15", "L_20",
      "E_15", "E_20",
      "C_Fin_15", "C_Fin_20",
      "Ka_15", "Ka_20",
      "Upp_15", "Upp_20",
      "VF_15", "VF_20",
      "NA_15", "NA_20",
      "NL_15", "NL_20"
    )
  )
  
  list(
    data = pcaData,
    percentVar = round(100 * attr(pcaData, "percentVar"))
  )
}

plot_pca <- function(pca, pc_x, pc_y, title) {
  ggplot(
    pca$data,
    aes(
      .data[[pc_x]],
      .data[[pc_y]],
      color = pop_temp,
      shape = region
    )
  ) +
    geom_point(size = 3) +
    scale_color_manual(
      values = pop_colors,
      name = "Population × temperature"
    ) +
    scale_shape_manual(
      name = "Region",
      values = c("North" = 16, "South" = 15, "East" = 14)
    ) +
    xlab(paste0(pc_x, ": ", pca$percentVar[1], "% variance")) +
    ylab(paste0(pc_y, ": ", pca$percentVar[2], "% variance")) +
    labs(title = title) +
    coord_fixed() +
    theme_minimal()
}

pca12 <- get_pca_data(vsd, c(1, 2))
plot_pca(
  pca12,
  "PC1", "PC2",
  "PCA of all samples: family, temperature, and region structure"
)

pca23 <- get_pca_data(vsd, c(2, 3))
plot_pca(
  pca23,
  "PC2", "PC3",
  "PCA of all samples: family, temperature, and region structure"
)

pca34 <- get_pca_data(vsd, c(3, 4))
plot_pca(
  pca34,
  "PC3", "PC4",
  "PCA of all samples: family, temperature, and region structure"
)

# ------------------------------------------------------------------------------
# Heatmap

get_top_DE_genes <- function(res, n = 400, padj_cutoff = 0.05, lfc_cutoff = 1) {
  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff &
               abs(res$log2FoldChange) >= lfc_cutoff, ]
  rownames(sig[order(sig$padj), ][1:min(n, nrow(sig)), ])
}

top_ES <- get_top_DE_genes(res_df_ES)
top_NS <- get_top_DE_genes(res_df_NS)
top_EN <- get_top_DE_genes(res_df_EN)

top_genes <- unique(c(top_ES, top_NS, top_EN))

# Column annotations (temperature is the focus; region/population are annotations only)
annotation_col <- as.data.frame(
  colData(dds)[, c("region", "temperature", "population")]
)

# Expression matrix for the already-defined top_genes
mat <- assay(vsd)[top_genes, ]

# Z-score per gene (pattern-focused for region)
mat_scaled <- t(scale(t(mat)))

# Annotation colors
ann_colors <- list(
  region = c(
    North = "#1f78b4",
    South = "#e31a1c",
    East = "#33a02c"
  ),
  temperature = c(
    "15" = "#a6cee3",
    "20" = "#fb9a99"
  ),
  population = c(
    "Ka"    = "#FFDE52",
    "Upp"   = "#FF7070",
    "NA"    = "#6f9fe2ff",
    "NL"    = "#59FFC0",
    "VF"    = "#6a3d9a",
    "C_Fin" = "#FF33E4",
    "E"     = "#B7FF5E",
    "L"     = "#00CDFF"
  )
)


ann_colors$row <- list(
  temp_effect = c(
    "Higher at 20" = "#b2182b",  # dark red
    "Higher at 15" = "#2166ac"   # dark blue
  )
)

# Heatmap: region-driven structure via correlation clustering
pheatmap(
  mat_scaled,
  main="Region-driven transcriptional structure",
  scale = "none",
  cluster_rows = T,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  fontsize_col = 8
)

# ------------------------------------------------------------------------------
# Heatmap
# EAST VS SOUTH

# Select significant genes
sig <- res_df_ES[!is.na(res_df_ES$padj) & res_df_ES$padj < 0.05, ]

# Top 10 upregulated genes
top_up <- sig[sig$log2FoldChange > 0, ]
top_up <- top_up[order(top_up$padj), ][1:400, ]

# Top 10 downregulated genes
top_down <- sig[sig$log2FoldChange < 0, ]
top_down <- top_down[order(top_down$padj), ][1:400, ]

# Combine genes
top_genes <- rbind(top_up, top_down)

# Column annotations (temperature is the focus; region/population are annotations only)
annotation_col <- as.data.frame(
  colData(dds)[, c("temperature", "region", "population")]
)

row_annotation <- data.frame(
  temp_effect = ifelse(
    res_df_ES[rownames(top_genes), "log2FoldChange"] > 0,
    "Higher at 20",
    "Higher at 15"
  )
)
rownames(row_annotation) <- rownames(top_genes)

# Expression matrix for the already-defined top_genes
mat <- assay(vsd)[rownames(top_genes), ]

# Z-score per gene (pattern-focused for temperature)
mat_scaled <- t(scale(t(mat)))

# Annotation colors
ann_colors <- list(
  region = c(
    East = "#33a02c",
    South = "#e31a1c"
  ),
  temperature = c(
    "15" = "#a6cee3",
    "20" = "#fb9a99"
  ),
  population = c(
    "Ka"    = "#FFDE52",
    "Upp"   = "#FF7070",
    "C_Fin" = "#FF33E4",
    "E"     = "#B7FF5E",
    "L"     = "#00CDFF"
  )
)

ann_colors$row <- list(
  temp_effect = c(
    "Higher at 20" = "#b2182b",  # dark red
    "Higher at 15" = "#2166ac"   # dark blue
  )
)

keep <- annotation_col$region %in% c("East", "South")

mat_scaled <- mat_scaled[, keep]
annotation_col <- annotation_col[keep, ]

# Heatmap: temperature-driven structure via correlation clustering
pheatmap(
  mat_scaled,
  main="Region-driven transcriptional structure (East vs South)",
  scale = "none",
  cluster_rows = T,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  fontsize_col = 8
)

# ------------------------------------------------------------------------------
# Heatmap NORTH VS SOUTH

# Select significant genes
sig <- res_df_NS[!is.na(res_df_NS$padj) & res_df_NS$padj < 0.05, ]

# Top 10 upregulated genes
top_up <- sig[sig$log2FoldChange > 0, ]
top_up <- top_up[order(top_up$padj), ][1:400, ]

# Top 10 downregulated genes
top_down <- sig[sig$log2FoldChange < 0, ]
top_down <- top_down[order(top_down$padj), ][1:400, ]

# Combine genes
top_genes <- rbind(top_up, top_down)

# Column annotations (temperature is the focus; region/population are annotations only)
annotation_col <- as.data.frame(
  colData(dds)[, c("temperature", "region", "population")]
)

row_annotation <- data.frame(
  temp_effect = ifelse(
    res_df_NS[rownames(top_genes), "log2FoldChange"] > 0,
    "Higher at 20",
    "Higher at 15"
  )
)
rownames(row_annotation) <- rownames(top_genes)

# Expression matrix for the already-defined top_genes
mat <- assay(vsd)[rownames(top_genes), ]

# Z-score per gene (pattern-focused for temperature)
mat_scaled <- t(scale(t(mat)))

# Annotation colors
ann_colors <- list(
  region = c(
    North = "#1f78b4",
    South = "#e31a1c"
  ),
  temperature = c(
    "15" = "#a6cee3",
    "20" = "#fb9a99"
  ),
  population = c(
    "NA"    = "#6f9fe2ff",
    "NL"    = "#59FFC0",
    "VF"    = "#6a3d9a",
    "Ka"    = "#FFDE52",
    "Upp"   = "#FF7070"
  )
)

ann_colors$row <- list(
  temp_effect = c(
    "Higher at 20" = "#b2182b",  # dark red
    "Higher at 15" = "#2166ac"   # dark blue
  )
)

keep <- annotation_col$region %in% c("North", "South")

mat_scaled <- mat_scaled[, keep]
annotation_col <- annotation_col[keep, ]

# Heatmap: temperature-driven structure via correlation clustering
pheatmap(
  mat_scaled,
  main="Region-driven transcriptional structure (North vs South)",
  scale = "none",
  cluster_rows = T,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  fontsize_col = 8
)

# ------------------------------------------------------------------------------
# Heatmap EAST VS NORTH

# Select significant genes
sig <- res_df_EN[!is.na(res_df_EN$padj) & res_df_EN$padj < 0.05, ]

# Top 10 upregulated genes
top_up <- sig[sig$log2FoldChange > 0, ]
top_up <- top_up[order(top_up$padj), ][1:400, ]

# Top 10 downregulated genes
top_down <- sig[sig$log2FoldChange < 0, ]
top_down <- top_down[order(top_down$padj), ][1:400, ]

# Combine genes
top_genes <- rbind(top_up, top_down)

# Column annotations (temperature is the focus; region/population are annotations only)
annotation_col <- as.data.frame(
  colData(dds)[, c("temperature", "region", "population")]
)

row_annotation <- data.frame(
  temp_effect = ifelse(
    res_df_EN[rownames(top_genes), "log2FoldChange"] > 0,
    "Higher at 20",
    "Higher at 15"
  )
)
rownames(row_annotation) <- rownames(top_genes)

# Expression matrix for the already-defined top_genes
mat <- assay(vsd)[rownames(top_genes), ]

# Z-score per gene (pattern-focused for temperature)
mat_scaled <- t(scale(t(mat)))

# Annotation colors
ann_colors <- list(
  region = c(
    East = "#33a02c",
    North = "#1f78b4"
  ),
  temperature = c(
    "15" = "#a6cee3",
    "20" = "#fb9a99"
  ),
  population = c(
    "NA"    = "#6f9fe2ff",
    "NL"    = "#59FFC0",
    "VF"    = "#6a3d9a",
    "C_Fin" = "#FF33E4",
    "E"     = "#B7FF5E",
    "L"     = "#00CDFF"
  )
)

ann_colors$row <- list(
  temp_effect = c(
    "Higher at 20" = "#b2182b",  # dark red
    "Higher at 15" = "#2166ac"   # dark blue
  )
)

keep <- annotation_col$region %in% c("East", "North")

mat_scaled <- mat_scaled[, keep]
annotation_col <- annotation_col[keep, ]

# Heatmap: temperature-driven structure via correlation clustering
pheatmap(
  mat_scaled,
  main="Region-driven transcriptional structure (East vs North)",
  scale = "none",
  cluster_rows = T,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  fontsize_col = 8
)



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

genes_ES <- rownames(res_df_ES)[
  !is.na(res_df_ES$padj) &
  res_df_ES$padj < 0.05 &
  abs(res_df_ES$log2FoldChange) >= 1
]

genes_NS <- rownames(res_df_NS)[
  !is.na(res_df_NS$padj) &
  res_df_NS$padj < 0.05 &
  abs(res_df_NS$log2FoldChange) >= 1
]

genes_EN <- rownames(res_df_EN)[
  !is.na(res_df_EN$padj) &
  res_df_EN$padj < 0.05 &
  abs(res_df_EN$log2FoldChange) >= 1
]

length(genes_ES)
length(genes_NS)
length(genes_EN)

length(intersect(genes_ES, genes_NS))
length(intersect(genes_ES, genes_EN))
length(intersect(genes_NS, genes_EN))

length(Reduce(intersect, list(genes_ES, genes_NS, genes_EN)))


library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    "East vs South" = genes_ES,
    "North vs South" = genes_NS,
    "East vs North" = genes_EN
  ),
  filename = NULL,
  fill = c("#e31a1c", "#1f78b4", "#33a02c"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2
)

grid.draw(venn.plot)


