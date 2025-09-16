import pandas as pd
from pyscenic.binarization import binarize
auc_mtx = pd.read_csv('dcm_hcm_healthy_auc.csv', index_col=0)
bin_mtx, thresholds = binarize(auc_mtx)
bin_mtx.to_csv('dcm_hcm_healthy_bin.csv')
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv('dcm_hcm_healthy_threshold.csv')
