import scanpy as sc
import pandas as pd
from gssnng import score_cells
import os




annda=sc.read_10x_mtx('./vcm_sample2')
meta=pd.read_csv('./vcm_sample2/meta.csv')
annda.obs=meta
sc.pp.normalize_total(annda, target_sum=1e4)
sc.pp.log1p(annda)


initial=len(annda.obs.columns)


score_cells.with_gene_sets(adata=annda,                            
                           gene_set_file='c2.cp.kegg_legacy.v2024.1.Hs.symbols_up.gmt', 
                           groupby='groups',                 
                           smooth_mode='connectivity',        
                           recompute_neighbors=32,             
                           score_method='singscore',          
                           method_params={'normalization':'theoretical'},  
                           ranked=True,                       
                           cores=10)
annda.obs.to_csv('./python_term_score_res/kegg_score_res_legacy.csv',index_label=False)

annda.obs=annda.obs.iloc[:,:initial]
					

total=os.listdir('/data/lixingyu/heart_origin_analysis/group_analysis2/dcm/vcm/go_bp')
					
for i in total:					   
	score_cells.with_gene_sets(adata=annda,                            
							   gene_set_file=f'/data/lixingyu/heart_origin_analysis/group_analysis2/dcm/vcm/go_bp/{i}', 
							   groupby='groups',                 
							   smooth_mode='connectivity',        
							   recompute_neighbors=32,             
							   score_method='singscore',          
							   method_params={'normalization':'theoretical'},  
							   ranked=True,                       
							   cores=10)
	cur_name=i.split('_')[1].split('.')[0]
	annda.obs.to_csv(f'./python_term_score_res/gobp_score_{cur_name}.csv',index_label=False)
	annda.obs=annda.obs.iloc[:,:initial]
