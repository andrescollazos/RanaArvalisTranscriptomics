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
  design = ~ region + temperature
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
save(list=ls(all=TRUE), file="06_Temp_Region.RData")
print("Saving Results: DONE")

print("-----------------------------------------------------------------")
warnings()

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
    title = "MA plot: temperature effect across all populations",
    x = "Mean of normalized counts",
    y = "Log2 fold change (shrunken)"
  ) +
  coord_cartesian(ylim = c(-6, 6)) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Volcano Plot

EnhancedVolcano(
  resLFC,
  lab = NA,  # no labels by default
  x = 'log2FoldChange',
  y = 'padj',
  
  xlim = c(-5, 5),
  ylim = c(0, 16),
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  pointSize = 1.8,
  labSize = 3.0,
  
  title = 'Volcano plot: temperature effect across all populations',
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
