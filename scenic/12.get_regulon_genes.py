import pickle
import pandas as pd
                    
with open('./dcm_hcm_healthy.regulons.dat','rb') as f:
    regulons=pickle.load(f)
total_tf=pd.read_csv('./intersect_hcm_compare_healthy.csv')
total_tf=total_tf.index
for i in total_tf:
    for j in regulons:
        if j.name==i:
            with open(f'./{i}_tf_hcm_healthy.txt','w') as f1:
                for k in list(j.gene2weight.keys()):
                    f1.write(k+'\n')




total_tf=pd.read_csv('./intersect_dcm_compare_healthy.csv')
total_tf=total_tf.index
for i in total_tf:
    for j in regulons:
        if j.name==i:
            with open(f'./{i}_tf_dcm_healthy.txt','w') as f1:
                for k in list(j.gene2weight.keys()):
                    f1.write(k+'\n')


total_tf=pd.read_csv('./intersect_dcm_compare_hcm.csv',index_col=0)
total_tf=total_tf.index
for i in total_tf:
    for j in regulons:
        if j.name==i:
            with open(f'./{i}_tf_dcm_hcm.txt','w') as f1:
                for k in list(j.gene2weight.keys()):
                    f1.write(k+'\n')
