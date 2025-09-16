import scanpy as sc
import pandas as pd
annda=sc.read_h5ad('/data/lixingyu/all_data_raw_with_cell_type_full.h5ad')
annda=annda[annda.obs.cell_type1=='vcm']
annda=annda[annda.obs.condition.isin(['DCM','HCM'])]
annda=annda[annda.obs.region.isin(['LV','IVS'])]


sample_num=annda.obs.loc[:,'sample'].value_counts()
sample_num=sample_num[sample_num>100]
for i in sample_num.index:
    cur_annda=annda[annda.obs.loc[:,'sample']==i]
    cur_raw_counts = cur_annda.X.toarray().T
    cur_raw_counts = pd.DataFrame(cur_raw_counts, index=cur_annda.var_names, columns=cur_annda.obs_names)
    cur_raw_counts.to_csv(f'/data/lixingyu/condition_analysis/dcm_hcm/vcm/{i}_rawcounts.csv')
