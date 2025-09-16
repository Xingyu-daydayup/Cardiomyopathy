library(glmmTMB)
library(bbmle)
number<-2000
success<-0
success_genes<-c()
lose_genes<-c()
pvalue <-c()
lose<-0
flag<-F
final_data<-read.csv('./model_input.csv',header = T)
total_genes<-colnames(final_data)[10:17215]
#total_genes<-colnames(my_data)[12:51]
res_coeff<-data.frame()
res_pval<-data.frame()
compute_deviance<-function(model,observed){
  predicted<-fitted(model)
  resid<-observed-predicted
  Poisson.dev<-function(y, mu){
    2*(y*log(ifelse(y == 0, 1, y/mu))-(y-mu))
  }
  residuals.deviance<- sqrt(Poisson.dev(observed,predicted))*ifelse(observed>predicted,1,-1)
  sum(residuals.deviance^2)
}
compute_p<-function(model,observed){
  1-pchisq(compute_deviance(model,observed), df.residual(model))
}

#final_data$condition<-factor(final_data$condition,levels = c('healthy','HCM'))
total_genes<-total_genes[(number-2000+1):number]
for(i in total_genes){
  observed <-final_data[,i]
  tryCatch({
    m1 <-glmmTMB(eval(parse(text=i))~condition+sex+age+race+offset(log(offset))+(1|donor),data = final_data,family = nbinom2)
  }, warning = function(w){
    if(flag==F){
      lose <<-lose+1
      pvalue<<-append(pvalue,'warn')
      lose_genes<<-append(lose_genes,i)
      flag<<-T}
  }, error = function(e){
    if(flag == F){
      lose <<-lose+1
      pvalue<<-append(pvalue,'error')
      lose_genes<<-append(lose_genes,i)
      flag<<-T}
  },finally = {
    if(flag==F){
      success <<-success+1
      success_genes<<-append(success_genes,i)
      p<-compute_p(m1,observed)
      pvalue <<-append(pvalue,p)
      p_val<-summary(m1)[[6]]$cond[,4]
      res_pval[rownames(summary(m1)[[6]]$cond),i]<-p_val
      res_coeff[rownames(summary(m1)[[6]]$cond),i]<-summary(m1)[[6]]$cond[,1]
      flag<<-T}
  })
  flag<<-F
}


write.table(data.frame(res_coeff),paste0('./model_res/',number,'_coeff.txt'))
write.table(data.frame(res_pval),paste0('./model_res/',number,'_pval.txt'))

write.csv(pvalue,paste0('./model_res/',number,'_pvals_model.txt'))
write.csv(success_genes,paste0('./model_res/',number,'_success_genes.txt'))
write.csv(lose_genes,paste0('./model_res/',number,'_lose_genes.txt'))
print(success)
print(lose)

