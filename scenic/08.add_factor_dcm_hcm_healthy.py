import pandas as pd

data=pd.read_csv('./dcm_hcm_healthy_bin.csv',index_col=0)
meta=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
meta=meta.loc[data.index,:]
data.loc[:,'condition']=meta.condition.values
data.loc[:,'region']=meta.region.values
data.loc[:,'sex']=meta.sex.values
data.loc[:,'age']=meta.age.values
data.loc[:,'race']=meta.race.values
data.to_csv('./dcm_hcm_healthy_bin_with_factor.csv',index_label=False)
