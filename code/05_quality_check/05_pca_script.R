args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript pca.R <expression_matrix_file> <samples_config_file> <output_pdf_name>")
}
matrix_file <- args[1]
samples_config_file <- args[2]
pdf_name <- paste0(args[3], ".pca.pdf")

library(ggplot2)
library(gridExtra)

## Load centered expression matrix (genes x samples)
print("Loading centered expression matrix")
data <- read.table(
  matrix_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)
print("Loading centered expression matrix: DONE")

## Load sample metadata
samples_data <- read.table(
  samples_config_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(samples_data) <- c(
  "sample_name",
  "replicate_name",
  "population",
  "temperature",
  "baltic_side"
)

## Ensure metadata order matches PCA columns
samples_data <- samples_data[match(colnames(data), samples_data$replicate_name), ]

## Run PCA (samples x genes)
print("Running PCA analysis")
pca <- prcomp(t(data), center = FALSE, scale. = FALSE)
## Percent variance
pc_pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100
print("Running PCA analysis: DONE")

print("Plotting PC's")
## Base dataframe (samples × PCs)
pca_df <- as.data.frame(pca$x)
pca_df$population   <- samples_data$population
pca_df$temperature  <- samples_data$temperature
pca_df$baltic_side  <- samples_data$baltic_side

## EXACT SAME COLOR LOGIC (verbatim)
pca_df$color <- ifelse(
  pca_df$baltic_side == "Swedish" & pca_df$temperature == 15, "violetred1",
  ifelse(
    pca_df$baltic_side == "Swedish" & pca_df$temperature == 20, "red4",
    ifelse(
      pca_df$baltic_side == "Non-Swedish" & pca_df$temperature == 15, "deepskyblue3",
      "blue"
    )
  )
)

## EXACT SAME SHAPE LOGIC
populations <- unique(samples_data$population)
pch_map <- setNames(seq_along(populations), populations)
pca_df$shape <- pch_map[pca_df$population]

plots <- list()
## GGPlot PCA loop
for (i in 1:5) {
  p <- ggplot(
    pca_df,
    aes(
      x = .data[[paste0("PC", i)]],
      y = .data[[paste0("PC", i + 1)]]
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
  geom_point(
    aes(color = color, shape = population),
    size = 3
  ) +

  ## COLOR LEGEND (using your existing colors)
  scale_color_identity(
    guide = "legend",
    name = "Baltic side x temperature",
    breaks = c("violetred1", "red4", "deepskyblue3", "blue"),
    labels = c(
      "Swedish 15°C",
      "Swedish 20°C",
      "Non-Swedish 15°C",
      "Non-Swedish 20°C"
    )
  ) +

  ## SHAPE LEGEND (unchanged)
  scale_shape_manual(
    values = pch_map,
    name = "Population",
    na.translate = FALSE
  ) +

  labs(
    x = paste0("PC", i, " (", round(pc_pct[i], 2), "%)"),
    y = paste0("PC", i + 1, " (", round(pc_pct[i + 1], 2), "%)")
  ) +

  ## FULL BORDER
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(colour = "black"),
    legend.position = "right"
  )

  plots[[i]] <- p
}

print("Saving PDF file")
ggsave(
  filename = pdf_name,
  plot = marrangeGrob(
    grobs = plots,
    nrow = 2,
    ncol = 1
  ),
  width = 10,
  height = 12
)
print("Saving PDF file: DONE")
