import pandas as pd
meta=pd.read_csv('/data/lixingyu/final_meta_for_model.csv')
meta=meta[meta.cell_type1=='vcm']
meta=meta[meta.condition.isin(['DCM','HCM'])]
meta=meta[meta.region.isin(['LV','IVS'])]



sample_num=meta.loc[:,'sample'].value_counts()
sample_num=sample_num[sample_num>100]
for i in sample_num.index:
    cur_meta=meta[meta.loc[:,'sample']==i]
    cur_meta=cur_meta.iloc[:,[-2]]
    cur_meta.to_csv(f'/data/lixingyu/condition_analysis/dcm_hcm/vcm/{i}_cluster_info.vcm.txt',header=False,index=True,index_label=False,sep='\t')
