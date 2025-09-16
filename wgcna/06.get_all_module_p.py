import pandas as pd
import os
total=os.listdir('./')

need=[]
for i in total:
	if i.endswith('_ans.csv'):
		need.append(i)

res=pd.read_csv(f'./{need[0]}',index_col=0)
res.columns=[need[0].split('_')[0]+'_pvalues']
for i in range(1,len(need)):
	cur_res=pd.read_csv(f'./{need[i]}',index_col=0)
	cur_res.columns=[need[i].split('_')[0]+'_pvalues']
	res=pd.concat([res,cur_res],axis=1)
res.to_csv('total_module_p.csv')
