import pandas as pd
import numpy as np
data=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
data=data[data.condition.isin(['DCM','healthy'])]
data=data[data.region.isin(['LV','IVS'])]

sample_num=data.loc[:,'sample'].value_counts()
sample_num=sample_num[sample_num>100]
data=data[data.loc[:,'sample'].isin(sample_num.index)]

data=data[data.loc[:,'sample']!='wrw_lv']

data=data.loc[:,['sample','condition','region','age','sex','data_name','cell_type1','race']]
data.to_csv('./meta_info.csv',index=False,index_label=False)
