library(ggplot2)
data<-read.csv('me_value_with_factor.csv')
sig_module<-read.csv('total_module_p.csv',row.names=1)
res<-colSums(sig_module<0.05)
res<-res[res>0]
for(i in names(res)){
  ggplot(data,aes(condition,eval(parse(text=strsplit(i,'_')[[1]][1]))))+
  geom_boxplot(aes(fill = condition))
  ggsave(paste0(strsplit(i,'_')[[1]][1],'.png'))
}
