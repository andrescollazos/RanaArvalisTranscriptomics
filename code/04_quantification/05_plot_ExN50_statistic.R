setwd('~/Documents/Studies/RanaArvalisTranscriptomics/results/04_quantification');

library(tidyverse)

# fixed filename
# filename <- "ExN50.gene.stats"
filename <- "ExN50.transcript.stats"

message(sprintf("parsing: %s", filename))

# read the data
alldata <- read.table(filename, header = TRUE, row.names = NULL)

# add sample name (optional but keeps the ggplot color consistent)
alldata$sample <- filename

# plot
p <- alldata %>%
  filter(Ex >= 30) %>%
  ggplot(aes(x = Ex, y = ExN50)) +
  geom_point(size = 2, color = "blue") +
  geom_line(color = "blue", linetype = "dashed") +
  xlim(c(30, 100)) +
  theme_bw()

print(p)
