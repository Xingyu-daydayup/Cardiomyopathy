import pandas as pd
res=pd.DataFrame()
for i in range(1,9):
    cur_data=pd.read_table(f'./model_res/{i*2000}_pval.txt',sep=' ')
    res=pd.concat([res,cur_data],axis=1)
res.to_csv('./model_res/total_factor_p.csv',index_label=False)

res=pd.DataFrame()
for i in range(1,9):
    cur_data = pd.read_table(f'./model_res/{i*2000}_coeff.txt', sep=' ')
    res = pd.concat([res, cur_data], axis=1)
res.to_csv('./model_res/total_coeff.csv',index_label=False)

