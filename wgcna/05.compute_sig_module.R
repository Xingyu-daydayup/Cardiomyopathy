library(MASS)
library(Matrix)
library(tidyverse)
library(multcomp)

data<-read.csv('./me_value_with_factor.csv')
colnames(data)

success_term_hcm<-c()
success_p_hcm<-c()

data$condition<-factor(data$condition,levels=c('healthy','HCM','DCM'))
for(i in 1:(length(colnames(data))-9)){
  if(colnames(data)[i] !='MEgrey'){
    model<-lm(data[,i]~age+race+sex+condition,data=data)
#    plot(model)
    res<-summary(glht(model, mcp(condition="Tukey")))
    ans<-data.frame(row.names=names(res$test$tstat))
    ans$pvalues<-res$test$pvalues
    write.csv(ans,paste0(colnames(data)[i],'_ans.csv'))
    
    
  }
  
}

