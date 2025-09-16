import pandas as pd

sig_tf = pd.read_csv('sig_tf_dcm_hcm_healthy.csv', index_col=0)
sig_tf.index = [i[:-5] for i in sig_tf.index]
coeff = pd.read_csv('tf_coeff_dcm_hcm_healthy.csv', index_col=0)
coeff.columns = [i[:-6] for i in coeff.columns]
coeff = coeff.T

# coeff>0是HCM高,
hcm_healthy = sig_tf[sig_tf.loc[:, 'HCM - healthy'] < 0.05]
cur_coeff = coeff.loc[hcm_healthy.index, :]
hcm_high = cur_coeff[cur_coeff.loc[:, 'HCM - healthy'] > 0]
healthy_high = cur_coeff[cur_coeff.loc[:, 'HCM - healthy'] < 0]
with open('./hcm_high_sig_tf_when_hcm_compare_healthy.txt','w') as f1:
    for i in hcm_high.index:
        f1.write(i+'\n')
with open('./healthy_high_sig_tf_when_hcm_compare_healthy.txt','w') as f1:
    for i in healthy_high.index:
        f1.write(i+'\n')

# coeff>0是DCM高
dcm_healthy = sig_tf[sig_tf.loc[:, 'DCM - healthy'] < 0.05]
cur_coeff = coeff.loc[dcm_healthy.index, :]
dcm_high = cur_coeff[cur_coeff.loc[:, 'DCM - healthy'] > 0]
healthy_high = cur_coeff[cur_coeff.loc[:, 'DCM - healthy'] < 0]
with open('./dcm_high_sig_tf_when_dcm_compare_healthy.txt','w') as f1:
    for i in dcm_high.index:
        f1.write(i + '\n')
with open('./healthy_high_sig_tf_when_dcm_compare_healthy.txt','w') as f1:
    for i in healthy_high.index:
        f1.write(i + '\n')

# coeff>0是DCM高
dcm_hcm = sig_tf[sig_tf.loc[:, 'DCM - HCM'] < 0.05]
cur_coeff = coeff.loc[dcm_hcm.index, :]
dcm_high = cur_coeff[cur_coeff.loc[:, 'DCM - HCM'] > 0]
hcm_high = cur_coeff[cur_coeff.loc[:, 'DCM - HCM'] < 0]
with open('./dcm_high_sig_tf_when_dcm_compare_hcm.txt','w') as f1:
    for i in dcm_high.index:
        f1.write(i + '\n')
with open('./hcm_high_sig_tf_when_dcm_compare_hcm.txt','w') as f1:
    for i in hcm_high.index:
        f1.write(i + '\n')

