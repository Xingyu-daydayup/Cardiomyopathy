library(tidyverse) 
library(magrittr)
library(WGCNA)
library(ggplot2)
data<-read.csv('wgcna_input.csv',check.names=F)
allowWGCNAThreads(4)
png('softpower.png')
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  data,             # <= Input data
 # blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
dev.off()
