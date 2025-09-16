import pandas as pd
auc_mtx=pd.read_csv('dcm_hcm_healthy_auc.csv',index_col=0)
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
meta=pd.read_csv('dcm_hcm_healthy_bin_with_factor.csv',index_col=0)
auc_mtx_Z.loc[:,'condition']=meta.condition
res=auc_mtx_Z.groupby('condition').agg('mean')
need=pd.read_csv('intersect_hcm_compare_healthy.csv')
need=res.loc[['HCM','healthy'],need.index]
need.to_csv('./hcm_healthy.csv',index_label=False)

need=pd.read_csv('intersect_dcm_compare_healthy.csv')
need=res.loc[['DCM','healthy'],need.index]
need.to_csv('./dcm_healthy.csv',index_label=False)

need=pd.read_csv('intersect_dcm_compare_hcm.csv')
need=res.loc[['DCM','HCM'],need.index]
need.to_csv('./dcm_hcm.csv',index_label=False)
