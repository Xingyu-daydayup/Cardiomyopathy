import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc
import seaborn as sns
import scanpy.external as sce
from anndata import  AnnData
import anndata as ad

csx_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/csx_lv/output_filtered.h5')
csx_lv.obs.index = ['csx_lv-'+i for i in csx_lv.obs.index]
csx_lv.obs['sample'] = 'csx_lv'
csx_lv.obs['condition'] = 'DCM'
csx_lv.obs['region'] = 'LV'
csx_lv.obs['cell_nuclei'] = 'nuclei'


csx_rv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/csx_rv/output_filtered.h5')
csx_rv.obs.index = ['csx_rv-'+i for i in csx_rv.obs.index]
csx_rv.obs['sample'] = 'csx_rv'
csx_rv.obs['condition'] = 'DCM'
csx_rv.obs['region'] = 'RV'
csx_rv.obs['cell_nuclei'] = 'nuclei'


lyj_avn = sc.read_10x_h5('/Lustre02/lixingyu/test_data/lyj_avn/output_filtered.h5')
lyj_avn.obs.index = ['lyj_avn-'+i for i in lyj_avn.obs.index]
lyj_avn.obs['sample'] = 'lyj_avn'
lyj_avn.obs['condition'] = 'ARVC'
lyj_avn.obs['region'] = 'AVN'
lyj_avn.obs['cell_nuclei'] = 'nuclei'


lyj_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/lyj_lv/output_filtered.h5')
lyj_lv.obs.index = ['lyj_lv-'+i for i in lyj_lv.obs.index]
lyj_lv.obs['sample'] = 'lyj_lv'
lyj_lv.obs['condition'] = 'ARVC'
lyj_lv.obs['region'] = 'LV'
lyj_lv.obs['cell_nuclei'] ='nuclei'


lyj_rv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/lyj_rv/output_filtered.h5')
lyj_rv.obs.index = ['lyj_rv-'+i for i in lyj_rv.obs.index]
lyj_rv.obs['sample'] = 'lyj_rv'
lyj_rv.obs['condition'] = 'ARVC'
lyj_rv.obs['region'] = 'RV'
lyj_rv.obs['cell_nuclei'] = 'nuclei'


rgf_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/rgf_lv/output_filtered.h5')
rgf_lv.obs.index = ['rgf_lv-'+i for i in rgf_lv.obs.index]
rgf_lv.obs['sample'] = 'rgf_lv'
rgf_lv.obs['condition'] = 'HCM'
rgf_lv.obs['region'] = 'LV'
rgf_lv.obs['cell_nuclei'] = 'nuclei'


scf_lv_sn = sc.read_10x_h5('/Lustre02/lixingyu/test_data/scf_lv_sn/output_filtered.h5')
scf_lv_sn.obs.index = ['scf_lv_sn-'+i for i in scf_lv_sn.obs.index]
scf_lv_sn.obs['sample'] = 'scf_lv_sn'
scf_lv_sn.obs['condition'] = 'DCM'
scf_lv_sn.obs['region'] = 'LV'
scf_lv_sn.obs['cell_nuclei'] = 'nuclei'


scf_ra = sc.read_10x_h5('/Lustre02/lixingyu/test_data/scf_ra/output_filtered.h5')
scf_ra.obs.index = ['scf_ra-'+i for i in scf_ra.obs.index]
scf_ra.obs['sample'] = 'scf_ra'
scf_ra.obs['condition'] = 'DCM'
scf_ra.obs['region'] = 'RA'
scf_ra.obs['cell_nuclei']='nuclei'


wqx_as = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wqx_as/output_filtered.h5')
wqx_as.obs.index = ['wqx_as-'+i for i in wqx_as.obs.index]
wqx_as.obs['sample'] = 'wqx_as'
wqx_as.obs['condition'] = 'HCM'
wqx_as.obs['region'] = 'AS'
wqx_as.obs['cell_nuclei'] = 'nuclei'


wqx_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wqx_lv/output_filtered.h5')
wqx_lv.obs.index = ['wqx_lv-'+i for i in wqx_lv.obs.index]
wqx_lv.obs['sample'] = 'wqx_lv'
wqx_lv.obs['condition'] = 'HCM'
wqx_lv.obs['region'] = 'LV'
wqx_lv.obs['cell_nuclei'] = 'nuclei'


wqx_rv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wqx_rv/output_filtered.h5')
wqx_rv.obs.index = ['wqx_rv-'+i for i in wqx_rv.obs.index]
wqx_rv.obs['sample'] = 'wqx_rv'
wqx_rv.obs['condition'] = 'HCM'
wqx_rv.obs['region'] = 'RV'
wqx_rv.obs['cell_nuclei'] = 'nuclei'


wrw_avn = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wrw_avn/output_filtered.h5')
wrw_avn.obs.index = ['wrw_avn-'+i for i in wrw_avn.obs.index]
wrw_avn.obs['sample'] = 'wrw_avn'
wrw_avn.obs['condition'] = 'DCM'
wrw_avn.obs['region'] = 'AVN'
wrw_avn.obs['cell_nuclei'] = 'nuclei'


wrw_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wrw_lv/output_filtered.h5')
wrw_lv.obs.index = ['wrw_lv-'+i for i in wrw_lv.obs.index]
wrw_lv.obs['sample'] = 'wrw_lv'
wrw_lv.obs['condition']='DCM'
wrw_lv.obs['region'] = 'LV'
wrw_lv.obs['cell_nuclei'] = 'nuclei'


wrw_rv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wrw_rv/output_filtered.h5')
wrw_rv.obs.index = ['wrw_rv-'+i for i in wrw_rv.obs.index]
wrw_rv.obs['sample'] = 'wrw_rv'
wrw_rv.obs['condition'] = 'DCM'
wrw_rv.obs['region'] = 'RV'
wrw_rv.obs['cell_nuclei'] = 'nuclei'


wrw_san = sc.read_10x_h5('/Lustre02/lixingyu/test_data/wrw_san/output_filtered.h5')
wrw_san.obs.index = ['wrw_san-'+i for i in wrw_san.obs.index]
wrw_san.obs['sample'] = 'wrw_san'
wrw_san.obs['condition'] = 'DCM'
wrw_san.obs['region'] = 'SAN'
wrw_san.obs['cell_nuclei'] = 'nuclei'


xby_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/xby_lv/output_filtered.h5')
xby_lv.obs.index = ['xby_lv-'+i for i in xby_lv.obs.index]
xby_lv.obs['sample'] = 'xby_lv'
xby_lv.obs['condition'] = 'DCM'
xby_lv.obs['region'] = 'LV'
xby_lv.obs['cell_nuclei'] = 'nuclei'


xby_ra = sc.read_10x_h5('/Lustre02/lixingyu/test_data/xby_ra/output_filtered.h5')
xby_ra.obs.index = ['xby_ra-'+i for i in xby_ra.obs.index]
xby_ra.obs['sample'] = 'xby_ra'
xby_ra.obs['condition'] = 'DCM'
xby_ra.obs['region'] = 'RA'
xby_ra.obs['cell_nuclei'] = 'nuclei'


hq_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/hq_lv/output_filtered.h5')
hq_lv.obs.index = ['hq_lv-'+i for i in hq_lv.obs.index]
hq_lv.obs['sample'] = 'hq_lv'
hq_lv.obs['condition'] = 'HCM'
hq_lv.obs['region'] = 'LV'
hq_lv.obs['cell_nuclei'] = 'nuclei'

lzb_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/lzb_lv/output_filtered.h5')
lzb_lv.obs.index = ['lzb_lv-'+i for i in lzb_lv.obs.index]
lzb_lv.obs['sample'] = 'lzb_lv'
lzb_lv.obs['condition'] = 'DCM'
lzb_lv.obs['region'] = 'LV'
lzb_lv.obs['cell_nuclei'] = 'nuclei'

jm_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/jm_lv/output_filtered.h5')
jm_lv.obs.index = ['jm_lv-'+i for i in jm_lv.obs.index]
jm_lv.obs['sample'] = 'jm_lv'
jm_lv.obs['condition'] = 'DCM'
jm_lv.obs['region'] = 'LV'
jm_lv.obs['cell_nuclei'] = 'nuclei'

zgw_alv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/zgw_alv/output_filtered.h5')
zgw_alv.obs.index = ['zgw_alv-'+i for i in zgw_alv.obs.index]
zgw_alv.obs['sample'] = 'zgw_alv'
zgw_alv.obs['condition'] = 'DCM'
zgw_alv.obs['region'] = 'ALV'
zgw_alv.obs['cell_nuclei'] = 'nuclei'

pqg_lv = sc.read_10x_h5('/Lustre02/lixingyu/test_data/pqg_lv/output_filtered.h5')
pqg_lv.obs.index = ['pqg_lv-'+i for i in pqg_lv.obs.index]
pqg_lv.obs['sample'] = 'pqg_lv'
pqg_lv.obs['condition'] = 'HCM'
pqg_lv.obs['region'] = 'LV'
pqg_lv.obs['cell_nuclei'] = 'nuclei'

ycz_lvot = sc.read_10x_h5('/Lustre02/lixingyu/test_data/ycz_lvot/output_filtered.h5')
ycz_lvot.obs.index = ['ycz_lvot-'+i for i in ycz_lvot.obs.index]
ycz_lvot.obs['sample'] = 'ycz_lvot'
ycz_lvot.obs['condition'] = 'HCM'
ycz_lvot.obs['region'] = 'LVOT'
ycz_lvot.obs['cell_nuclei'] = 'nuclei'

dbl_score_thresh=0.3
total_obj = [csx_rv,lyj_avn,lyj_lv,lyj_rv,rgf_lv,scf_lv_sn,scf_ra,wqx_as,wqx_lv,wqx_rv,wrw_avn,wrw_lv,wrw_rv,wrw_san,xby_lv,xby_ra,hq_lv,lzb_lv,jm_lv,zgw_alv,pqg_lv,ycz_lvot]
cur_obj=csx_lv
cur_obj.var_names_make_unique() 
scrub = scr.Scrublet(cur_obj.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    cur_obj.obs['scrublet_score'] = doublet_scores
    prop = np.sum(cur_obj.obs['scrublet_score']>=dbl_score_thresh)/len(cur_obj)
    print(f'proportion of doublets: {prop}')
    cur_obj = cur_obj[cur_obj.obs['scrublet_score']<dbl_score_thresh]


for i in total_obj:
    i.var_names_make_unique() 
    scrub = scr.Scrublet(i.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    i.obs['scrublet_score'] = doublet_scores
    prop = np.sum(i.obs['scrublet_score']>=dbl_score_thresh)/len(merge_data)
    print(f'proportion of doublets: {prop}')
    i = i[i.obs['scrublet_score']<dbl_score_thresh]
    cur_obj=ad.concat([cur_obj,i],join='outer')    
merge_data=cur_obj    
sc.pp.filter_cells(merge_data, min_genes=200,inplace=True)
sc.pp.filter_genes(merge_data, min_cells=3,inplace=True)
merge_data.var['mt'] = merge_data.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(merge_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
merge_data= merge_data[merge_data.obs.n_genes_by_counts < 10000, :]
merge_data = merge_data[merge_data.obs.pct_counts_mt < 5, :]
merge_data.write('qc.h5ad')
