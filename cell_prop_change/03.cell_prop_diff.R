library(speckle)
library(limma)
library(edgeR)
library(pheatmap)

cell_prop<-read.delim('cell_prop.txt',row.names = 1)
sample_num<-read.csv('sample_num.csv',header = T)
prop.list <- convertDataToList(cell_prop,data.type="proportions", transform="logit",
                               scale.fac=sample_num$sample)
sampleinfo <- read.csv("sampleinfo.csv")

designAS <- model.matrix(~0+sampleinfo$condition + sampleinfo$sex+sampleinfo$age+sampleinfo$race)

new_name<-c()
for(i in colnames(designAS)){
  cur_name<-strsplit(i, "$",fixed=T)[[1]][2]
  new_name<-append(new_name,cur_name)
}

colnames(designAS) <- new_name



#colnames(designAS) <- c("DCM","healthy","sex",'age','IVS','LA','LV','RA','RV','white','yellow')
mycontr <- makeContrasts(conditionDCM-conditionhealthy, levels=designAS)
res<-propeller.ttest(prop.list = prop.list,design = designAS, contrasts = mycontr,robust=TRUE,trend=FALSE,sort=TRUE)
write.table(res,'cell_prop_res.txt',quote = F)
