library(CellChat)
library(patchwork)
hcm <- readRDS('hcm_res.rds')
hcm <- netAnalysis_computeCentrality(hcm, slot.name = "netP") 
dcm<-readRDS('dcm_res.rds')
dcm <- netAnalysis_computeCentrality(dcm, slot.name = "netP") 
object.list <- list(HCM = hcm, DCM = dcm)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg3<-gg1 + gg2
ggsave('./figures/dcm_hcm/interactions_compare.png')


pdf('./figures/dcm_hcm/num_change.pdf')
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
dev.off()
pdf('./figures/dcm_hcm/weight_change.pdf')
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()




cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = c('uwot'))
cellchat <- netClustering(cellchat, type = "functional")
pdf('./figures/dcm_hcm/signaling_similarity_functional.pdf')
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()


cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = c('uwot'))
cellchat <- netClustering(cellchat, type = "structural")
pdf('./figures/dcm_hcm/signaling_similarity_structural.pdf')
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()



gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg3<-gg1 + gg2
ggsave('./figures/dcm_hcm/bar_plot.pdf')


library(ComplexHeatmap)


pdf('./figures/dcm_hcm/outgoing_signal.pdf',height=10)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


pdf('./figures/dcm_hcm/incoming_signal.pdf',height=10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf('./figures/dcm_hcm/overall_signal.pdf',height=10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
