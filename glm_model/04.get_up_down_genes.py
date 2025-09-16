import pandas as pd
coeff=pd.read_csv('./model_res/total_coeff.csv')
sig_genes=pd.read_table('./model_res/condition_genes_sig.txt',header=None)
sig_genes=sig_genes.iloc[:,0].values

coeff=coeff.loc['conditionHCM',:]
coeff=coeff[sig_genes]
up_genes=coeff[coeff>0]
down_genes=coeff[coeff<0]

with open('./model_res/up_genes.txt','w') as f1:
    for i in up_genes.index:
        f1.write(i+'\n')
with open('./model_res/down_genes.txt','w') as f1:
    for i in down_genes.index:
        f1.write(i+'\n')
