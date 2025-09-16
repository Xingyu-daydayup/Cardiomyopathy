library(CellChat)
load('./cell_chat_work_dcm.RData')
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)
groupSize <- as.numeric(table(cellChat@idents))
pathways.show <- cellChat@netP$pathways
for(i in pathways.show){
  pdf(paste0('./figures/dcm/',i,'_signaling_network.pdf'))
  par(mfrow=c(1,1),xpd=TRUE)
  netVisual_chord_cell(cellChat, signaling = i, title.name = paste0(i, " signaling network"))
  dev.off() 
}
for(i in pathways.show){
pdf(paste0('./figures/dcm/',i,'_signaling_network_circle.pdf'))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_aggregate(cellChat, signaling =i, layout = "circle")
dev.off()
}




