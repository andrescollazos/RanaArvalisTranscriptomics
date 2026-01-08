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

# Filter for Swedish-only samples
col_swe <- col_data[col_data$baltic_side == "Swe", ]
col_swe <- droplevels(col_swe)
# Subset count matrix to the same samples
cts <- counts_data[, rownames(col_swe)]
# Round matrix
cts <- round(cts)
# Check matrix and coldata structure
cat("Check matrix and coldata structure. Is it correct?", ncol(cts) == nrow(col_swe), "\n")

# Construct DESeq2 dataset
# Put the variable of interest at the end of the formula and make sure the control level is the first level.
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = col_swe,
  design = ~ population + family + temperature
)

# Pre filtering
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# Set Factor Level
dds$temperature <- relevel(dds$temperature, ref = "15")

# Run DESeq
print("Running DESeq")
dds <- DESeq(dds)
print("Obtaining results")
res <- results(dds, name = "temperature_20_vs_15")
print("Running DESeq: DONE")

print("Calculating Shrinkage of effect size")
resLFC <- lfcShrink(dds, coef = "temperature_20_vs_15", type = "apeglm")

print("Saving Results")
save(list=ls(all=TRUE), file="01_swe.RData")
print("Saving Results: DONE")


summary(res)
summary(resLFC)

# Convert to data frame
res_df <- as.data.frame(resLFC)

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

ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = diffexp_sig)) +
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
    title = "MA plot: temperature effect in Swedish populations",
    x = "Mean of normalized counts",
    y = "Log2 fold change (shrunken)"
  ) +
  coord_cartesian(ylim = c(-6, 6)) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Volcano Plot

EnhancedVolcano(
  res,
  lab = NA,  # no labels by default
  x = 'log2FoldChange',
  y = 'padj',
  
  ylim = c(0, 20),
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  pointSize = 1.8,
  labSize = 3.0,
  
  title = 'Volcano plot: temperature effect in Swedish populations',
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
      "Ka_15", "Ka_20",
      "Upp_15", "Upp_20",
      "NA_15", "NA_20",
      "NL_15", "NL_20",
      "VF_15", "VF_20"
    )
  )
  
  list(
    data = pcaData,
    percentVar = round(100 * attr(pcaData, "percentVar"))
  )
}

plot_pca <- function(pca, pc_x, pc_y, xlab_idx, ylab_idx, title) {
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
      values = c("North" = 16, "South" = 15)
    ) +
    xlab(paste0(pc_x, ": ", pca$percentVar[xlab_idx], "% variance")) +
    ylab(paste0(pc_y, ": ", pca$percentVar[ylab_idx], "% variance")) +
    labs(title = title) +
    coord_fixed() +
    theme_minimal()
}

pca12 <- get_pca_data(vsd, c(1, 2))
plot_pca(
  pca12,
  "PC1", "PC2",
  1, 2,
  "PCA of Swedish samples: population, temperature, and region structure"
)

pca23 <- get_pca_data(vsd, c(2, 3))
plot_pca(
  pca23,
  "PC2", "PC3",
  1, 2,
  "PCA of Swedish samples: population, temperature, and region structure"
)
# ------------------------------------------------------------------------------
# Heatmap

# Select significant genes
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

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
    res_df[rownames(top_genes), "log2FoldChange"] > 0,
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
    Ka  = "#fdbf6f",
    Upp = "#ff7f00",
    "NA"  = "#b2df8a",
    NL  = "#33a02c",
    VF  = "#6a3d9a"
  )
)

ann_colors$row <- list(
  temp_effect = c(
    "Higher at 20" = "#b2182b",  # dark red
    "Higher at 15" = "#2166ac"   # dark blue
  )
)

# Heatmap: temperature-driven structure via correlation clustering
pheatmap(
  mat_scaled,
  main="Temperature-driven transcriptional structure",
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

# Expression matrix (absolute scale)
mat_mag <- assay(vsd)[rownames(top_genes), ]
pheatmap(
  mat_mag,
  main="Absolute expression",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  fontsize_col = 8
)

# Explicit column order: population → temperature
for (pop in unique(annotation_col$population)) {
  
  idx <- annotation_col$population == pop
  
  mat_pop <- mat_scaled[, idx]
  
  # Remove genes with zero variance within this population
  keep <- apply(mat_pop, 1, sd, na.rm = TRUE) > 0
  mat_pop <- mat_pop[keep, ]
  
  pheatmap(
    mat_pop,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    annotation_col = annotation_col[idx, ],
    annotation_colors = ann_colors,
    annotation_row = row_annotation,
    show_rownames = FALSE,
    fontsize_col = 8,
    main = paste("Temperature response in population", pop)
  )
}