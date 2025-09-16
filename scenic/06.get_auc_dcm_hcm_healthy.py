import pandas as pd
import scanpy as sc
import pickle

from pyscenic.aucell import aucell

with open('dcm_hcm_healthy.regulons.dat','rb') as f:
    regulons = pickle.load(f)
annda=sc.read_loom('DCM_HCM_healthy_data.loom')
exp_mtx=pd.DataFrame(annda.X.toarray(),index=annda.obs_names,columns=annda.var_names)
auc_mtx = aucell(exp_mtx, regulons, num_workers=20)
auc_mtx.to_csv('./dcm_hcm_healthy_auc.csv')
