import pandas as pd
me=pd.read_table('me_value.txt',sep=' ')
meta=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
meta=meta.loc[:,['sample','condition','region','age','race','donor','data_name','sex']]
meta=meta.drop_duplicates()
meta.set_index('sample',inplace=True)

right_order=[i[4:] for i in me.index]

meta=meta.loc[right_order,:]

me.loc[:,'sample']=meta.index
me.loc[:,'condition']=meta.loc[:,'condition'].values
me.loc[:,'region']=meta.loc[:,'region'].values
me.loc[:,'age']=meta.loc[:,'age'].values
me.loc[:,'race']=meta.loc[:,'race'].values
me.loc[:,'donor']=meta.loc[:,'donor'].values
me.loc[:,'data_name']=meta.loc[:,'data_name'].values
me.loc[:,'sex']=meta.loc[:,'sex'].values
me.to_csv('me_value_with_factor.csv',index_label=False)
