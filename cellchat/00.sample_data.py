import pandas as pd
import scanpy as sc
from scipy import io
meta=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
annda=sc.read_h5ad('/data/lixingyu/all_data_raw_with_cell_type_full.h5ad')
meta=meta[meta.condition.isin(['DCM','healthy','HCM'])]
meta=meta[meta.region.isin(['LV','IVS'])]
meta=meta[meta.cell_type1!='acm']
dcm=meta[meta.condition=='DCM']
healthy=meta[meta.condition=='healthy']
hcm=meta[meta.condition=='HCM']
need_dcm=[]
need_healthy=[]
need_hcm=[]
for i in meta.cell_type1.value_counts().index:
        cur_hcm=hcm[hcm.cell_type1==i]
        cur_dcm=dcm[dcm.cell_type1==i]
        cur_healthy=healthy[healthy.cell_type1==i]
        cur_hcm=cur_hcm.sample(frac=0.1,replace=False)
        cur_dcm=cur_dcm.sample(frac=0.1,replace=False)
        cur_healthy=cur_healthy.sample(frac=0.1,replace=False)
        need_dcm.extend(list(cur_dcm.index))
        need_healthy.extend(list(cur_healthy.index))
        need_hcm.extend(list(cur_hcm.index))
annda_dcm=annda[need_dcm,:]
annda_healthy=annda[need_healthy,:]
annda_hcm=annda[need_hcm,:]
with open('./dcm_sample/barcodes.tsv','w') as f1:
    for i in annda_dcm.obs_names:
        f1.write(i+'\n')
with open('./dcm_sample/features.tsv','w') as f1:
    for i in ['\t'.join([x,x,'Gene Expression']) for x in annda_dcm.var_names]:
        f1.write(i+'\n')
io.mmwrite('./dcm_sample/matrix.mtx',annda_dcm.X.T)
annda_dcm.obs.to_csv('./dcm_sample/meta.csv',index_label=False)




with open('./hcm_sample/barcodes.tsv','w') as f1:
    for i in annda_hcm.obs_names:
        f1.write(i+'\n')
with open('./hcm_sample/features.tsv','w') as f1:
    for i in ['\t'.join([x,x,'Gene Expression']) for x in annda_hcm.var_names]:
        f1.write(i+'\n')
io.mmwrite('./hcm_sample/matrix.mtx',annda_hcm.X.T)
annda_hcm.obs.to_csv('./hcm_sample/meta.csv',index_label=False)





with open('./healthy_sample/barcodes.tsv','w') as f1:
    for i in annda_healthy.obs_names:
        f1.write(i+'\n')
with open('./healthy_sample/features.tsv','w') as f1:
    for i in ['\t'.join([x,x,'Gene Expression']) for x in annda_healthy.var_names]:
        f1.write(i+'\n')
io.mmwrite('./healthy_sample/matrix.mtx',annda_healthy.X.T)
annda_healthy.obs.to_csv('./healthy_sample/meta.csv',index_label=False)
