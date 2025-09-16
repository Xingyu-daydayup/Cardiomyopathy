import scanpy as sc
import loompy as lp
import numpy as np
annda=sc.read_h5ad('/data/lixingyu/all_data_raw_with_cell_type_full.h5ad')
annda=annda[annda.obs.condition.isin(['DCM','HCM','healthy'])]
annda=annda[annda.obs.cell_type1=='vcm']
annda=annda[annda.obs.region.isin(['LV','IVS'])]
meta=annda.obs

sample_num=annda.obs.loc[:,'sample'].value_counts()
sample_num=sample_num[sample_num>100]

sample_index=[]
for cur_sample in sample_num.index:
    if sample_num[cur_sample] > 400:
        cur_meta=meta[meta.loc[:,'sample']==cur_sample]
        cur_meta=cur_meta.sample(400,replace=False)
        sample_index.extend(list(cur_meta.index))
    else:
        cur_meta = meta[meta.loc[:, 'sample'] == cur_sample]
        sample_index.extend(list(cur_meta.index))
annda=annda[sample_index,:]


row_attrs = {
            "Gene": np.array(annda.var.index) ,
            }
col_attrs = {
            "CellID":  np.array(annda.obs.index) ,
                "nGene": np.array( np.sum(annda.X.transpose()>0 , axis=0)).flatten() ,
                    "nUMI": np.array( np.sum(annda.X.transpose() , axis=0)).flatten() ,
                    }

lp.create('./DCM_HCM_healthy_data.loom', annda.X.transpose(), row_attrs, col_attrs )
