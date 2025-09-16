#library(ggpubr)
#data<-read.csv('C:/Users/DELL/Desktop/meta.csv')
#table(data$Group)
#my_comparisons <- list(c('CTL','KD'),c('CTL','KD_A'))
#ggplot(data, aes(x=Group, y=Neutrophil_Activation1,fill=Group)) +
#  geom_violin(trim=FALSE)+geom_boxplot(width=0.1)+stat_compare_means(comparisons = my_comparisons)+
#  stat_compare_means(label.y = 2.5)
#ggsave('test.pdf')


library(plyr)
data<-read.csv('./bulk_data_gobp_res.csv',row.names=1,comment.char = '*',check.names = F)

raw<-data$groups
raw1<-c('group1','group3','group4','group5','healthy')
new1<-c(1,2,3,4,5)
data$groups_continues<-mapvalues(raw,raw1,new1)


success_term<-c()
p_value<-c()
flag<-F
for(i in 2:(length(colnames(data))-5)){
  tryCatch({
    model1<-lm(data[,i]~data$age+data$race+data$sex+data$condition)
    model2<-lm(data[,i]~data$age+data$race+data$sex+data$condition+data$groups_continues)
  },error=function(e){
    if(flag==F){
      flag=T
    }
  },finally = {
    if(flag==F){
      res<-anova(model1, model2)
      if(res$`Pr(>F)`[2]<0.05){
        success_term<-append(success_term,colnames(data)[i])
        p_value<-append(p_value,res$`Pr(>F)`[2])
      }
      
    }
  })
  flag<-F
  
  
}


p_value<-p.adjust(p_value,method = 'BH')
res<-data.frame(row.names=success_term)
res$p_adj<-p_value
write.csv(res,'gobp_res_continues.csv')

