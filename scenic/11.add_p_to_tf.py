import pandas as pd
data=pd.read_csv('intersect_hcm_compare_healthy.csv')
meta=pd.read_csv('sig_tf_dcm_hcm_healthy.csv',index_col=0)
meta.index=[i[:-5]  for i in meta.index]
cur_meta=meta.loc[data.index,'HCM - healthy']
data.loc[:,'p_adj']=cur_meta
data.to_csv('./intersect_hcm_compare_healthy.csv',index_label=False)


data=pd.read_csv('intersect_dcm_compare_healthy.csv')
meta=pd.read_csv('sig_tf_dcm_hcm_healthy.csv',index_col=0)
meta.index=[i[:-5]  for i in meta.index]
cur_meta=meta.loc[data.index,'DCM - healthy']
data.loc[:,'p_adj']=cur_meta
data.to_csv('./intersect_dcm_compare_healthy.csv',index_label=False)

data=pd.read_csv('intersect_dcm_compare_hcm.csv')
meta=pd.read_csv('sig_tf_dcm_hcm_healthy.csv',index_col=0)
meta.index=[i[:-5]  for i in meta.index]
cur_meta=meta.loc[data.index,'DCM - HCM']
data.loc[:,'p_adj']=cur_meta
data.to_csv('./intersect_dcm_compare_hcm.csv',index_label=False)
