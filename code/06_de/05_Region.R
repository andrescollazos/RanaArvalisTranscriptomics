args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript diff_expression_swe.R <raw_count_matrix_file> <col_data_file>")
}
matrix_file <- args[1]
col_data_file <- args[2]

library("DESeq2")

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
  design = ~ temperature + population + family + region
)

# Pre filtering
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# Set Factor Level
dds$region <- relevel(dds$region, ref = "South")

# Run DESeq
print("Running DESeq")
dds <- DESeq(dds)
print("Obtaining results")
res_East_South <- results(dds, name = "region_East_vs_South")
res_North_South <- results(dds, name = "region_North_vs_South")
res_East_North <- results(dds, contrast = c("region", "East", "North"))
print("Running DESeq: DONE")

print("Calculating Shrinkage of effect size")
resLFC_East_South  <- lfcShrink(dds, coef = "region_East_vs_South",  type = "apeglm")
resLFC_North_South <- lfcShrink(dds, coef = "region_North_vs_South", type = "apeglm")
resLFC_East_North <- lfcShrink(dds, contrast = c("region", "East", "North"), type = "apeglm")

print("Saving Results")
save(list=ls(all=TRUE), file="05_Region.RData")
print("Saving Results: DONE")