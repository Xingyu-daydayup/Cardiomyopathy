library(stats)

library(stats)
library(MASS)
library(Matrix)
library(tidyverse)
library(multcomp)


data<-read.csv('./dcm_hcm_healthy_bin_with_factor.csv')
data$condition<-factor(data$condition,levels=c('healthy','HCM','DCM'))
model <- glm(data[,1] ~ condition + age+sex+race,data=data,family = binomial)
res<- summary(glht(model, mcp(condition="Tukey")))
ans<-data.frame(row.names=names(res$test$tstat))
ans[,paste0(colnames(data)[1],'_pval')]<-res$test$pvalues

coeff<-data.frame(row.names=names(res$test$tstat))
coeff[,paste0(colnames(data)[1],'_coeff')]<-res$test$coefficients



for(i in 2:(length(colnames(data))-5)){
  model <- glm(data[,i] ~ condition + age+sex+race,data=data,family = binomial)
  res<- summary(glht(model, mcp(condition="Tukey")))
  ans[,paste0(colnames(data)[i],'_pval')]<-res$test$pvalues
  coeff[,paste0(colnames(data)[i],'_coeff')]<-res$test$coefficients
  
  
}

ans<-t(ans)
ans[,1]<-p.adjust(ans[,1],method='BH')
ans[,2]<-p.adjust(ans[,2],method='BH')
ans[,3]<-p.adjust(ans[,3],method='BH')
write.csv(ans,'./sig_tf_dcm_hcm_healthy.csv')
write.csv(coeff,'./tf_coeff_dcm_hcm_healthy.csv')
