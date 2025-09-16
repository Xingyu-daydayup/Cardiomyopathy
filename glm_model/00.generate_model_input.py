import pandas as pd
import numpy as np
data=pd.read_table('../combined_pseudobulk.vcm.sample.filtered.PG.xls')
meta=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
meta.fillna(52,inplace=True)
meta=meta.loc[:,['sample','condition','age','race','donor','sex','region','data_name']]
meta=meta.drop_duplicates()
meta.set_index('sample',inplace=True)
sample_name=[]
for i in data.columns:
    sample_name.append('_'.join(i.split('_')[1:]))
meta=meta.loc[sample_name,:]

res=pd.DataFrame()
res.loc[:,'sample']=sample_name
res.loc[:,'condition']=meta.condition.values
res.loc[:,'age']=meta.age.values
res.loc[:,'race']=meta.race.values
res.loc[:,'donor']=meta.donor.values
res.loc[:,'sex']=meta.sex.values
res.loc[:,'region']=meta.region.values
res.loc[:,'data_name']=meta.data_name.values
res.loc[:,'offset']=np.sum(data).values

data=data.T

data.index=res.index
res=pd.concat([res,data],axis=1)
res.to_csv('./model_input.csv',index=False,index_label=False)
