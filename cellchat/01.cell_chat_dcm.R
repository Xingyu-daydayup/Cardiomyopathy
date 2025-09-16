library(Seurat)
library(CellChat)
library(patchwork)
options(future.globals.maxSize= 891289600)
pbmc.data <- Read10X(data.dir = "./dcm_sample")
meta<-read.csv('./dcm_sample/meta.csv')
data <- CreateSeuratObject(counts = pbmc.data,min.cells = 3, min.features = 200,meta=meta)
data <- NormalizeData(data)


cellChat <- createCellChat(object = data, group.by = "cell_type1", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)
cellChat@DB <- CellChatDB.use

cellChat <- subsetData(cellChat)
future::plan("multisession", workers = 10)
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)
df.net <- subsetCommunication(cellChat)
save.image('cell_chat_work_dcm.RData')
