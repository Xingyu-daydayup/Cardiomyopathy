import pandas as pd
import numpy as np
meta=pd.read_csv('./meta_info.csv')
samples_num=meta.loc[:,'sample'].value_counts()
samples_num.to_csv('./sample_num.csv',index_label=False)

cell_types=meta.cell_type1.value_counts().index
samples=meta.loc[:,'sample'].value_counts().index
res=pd.DataFrame(None,index=cell_types,columns=samples)


for i in samples:
    cur_sample_data=meta[meta.loc[:,'sample']==i]
    cur_cell_num=cur_sample_data.cell_type1.value_counts()
    cur_cell_num=cur_cell_num/np.sum(cur_cell_num)
    res.loc[cur_cell_num.index,i]=cur_cell_num.values
res.fillna(0,inplace=True)
res.to_csv('./cell_prop.txt',sep='\t',index_label=False)

meta_info=meta.loc[:,['sample','condition','region','sex','age','data_name','race']]
meta_info=meta_info.drop_duplicates()
meta_info=meta_info.reset_index()
meta_info.drop('index',axis=1,inplace=True)
meta_info.set_index('sample',inplace=True)
meta_info=meta_info.loc[res.columns,:]
meta_info.to_csv('sampleinfo.csv',index_label=False)
