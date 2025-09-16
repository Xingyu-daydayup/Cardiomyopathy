library(tidyverse) 
library(magrittr)
library(WGCNA)
library(ggplot2)
data<-read.csv('wgcna_input.csv',header=T,row.names=1,check.names=F,comment.char='*')
#data<-t(data)
picked_power = 6

temp_cor <- cor       
cor <- WGCNA::cor

netwk <- blockwiseModules(data,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = ncol(data),
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor
mergedColors = labels2colors(netwk$colors)
png(filename = 'dendrograms.png')
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

MEs0 <- moduleEigengenes(data, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$treatment = row.names(MEs0)
write.table(MEs0,'me_value.txt')
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
ggsave('model_sample_heatmap.png')
save.image("my_workspace.RData")
