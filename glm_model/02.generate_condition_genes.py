import pandas as pd
data=pd.read_csv('./model_res/total_factor_p.csv')
region=data.loc['conditionHCM',:]
region=region[~region.isna()]
region.to_csv('./model_res/condition_genes_p.csv',index_label=False)
