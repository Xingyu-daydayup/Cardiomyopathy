import pandas as pd
import numpy as np
data=pd.read_table('gene_modules.txt')

sig_module=pd.read_csv('./total_module_p.csv',index_col=0)
res=np.sum(sig_module<0.05)
res=res[res!=0]
for i in res.index:
    cur_name=i.split('_')[0][2:]
    cur_data=data[data.colors==cur_name]
    with open(f'./{cur_name}_genes.txt','w') as f1:
        for i in cur_data.gene_id.values:
            f1.write(i+'\n')
