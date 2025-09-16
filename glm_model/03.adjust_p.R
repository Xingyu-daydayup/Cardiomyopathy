data<-read.csv('./model_res/condition_genes_p.csv')
data$res<-p.adjust(data$conditionHCM,method='BH')
data<-data[data$res<0.05,]
write.table(rownames(data),'./model_res/condition_genes_sig.txt',quote = F,row.names = F,col.names = F)
