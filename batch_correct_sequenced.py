import scanpy as sc
import anndata as ad
import scanpy as sc
import seaborn as sns
import scanpy.external as sce
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from anndata import  AnnData
import os
import resource
soft_limit =1024* 1024 * 1024 * 1024  # 1GB*1024
hard_limit =1024* 1024 * 1024 * 1024 
resource.setrlimit(resource.RLIMIT_AS, (soft_limit, hard_limit))
os.environ['OPENBLAS_NUM_THREADS'] = '1'

my_nuclei_data = sc.read_h5ad('/Lustre02/lixingyu/qc.h5ad')


merge_data =my_nuclei_data


sc.pp.normalize_total(merge_data, target_sum=1e4)
sc.pp.log1p(merge_data)
sc.pp.highly_variable_genes(merge_data, min_mean=0.0125, max_mean=3, min_disp=0.25)
# This saves the original set of genes 
merge_data.raw = merge_data

merge_data = merge_data[:,merge_data.var.highly_variable]
sc.pp.regress_out(merge_data, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(merge_data, max_value=10)
sc.tl.pca(merge_data, svd_solver='arpack')
sc.pp.neighbors(merge_data, n_neighbors=10, n_pcs=50)

def one_col_lgd(umap):
    legend = umap.legend(bbox_to_anchor=[1.00, 0.5],
    loc='center left', ncol=1, prop={'size': 6})
    legend.get_frame().set_linewidth(0.0)
    for handle in legend.legendHandles:
        handle.set_sizes([25.0])
    return legend


sce.pp.harmony_integrate(merge_data, 'sample')
merge_data.obsm['X_pca'] = merge_data.obsm['X_pca_harmony']
sc.pp.neighbors(merge_data, n_neighbors=10, n_pcs=50)
sc.tl.umap(merge_data)
sc.tl.leiden(merge_data, resolution=0.2,key_added='leiden0.2')
sc.tl.leiden(merge_data, resolution=0.3,key_added='leiden0.3')
sc.tl.leiden(merge_data, resolution=0.4,key_added='leiden0.4')
sc.tl.leiden(merge_data, resolution=0.5,key_added='leiden0.5')
sc.tl.leiden(merge_data, resolution=0.6,key_added='leiden0.6')
sc.tl.leiden(merge_data, resolution=0.7,key_added='leiden0.7')

leiden_umap = sc.pl.umap(merge_data, color=['leiden0.4'],
show=False, palette=sns.color_palette("husl", 50),
    legend_fontsize=6, frameon=True, title='Leiden')

lgd = one_col_lgd(leiden_umap)

fig = leiden_umap.get_figure()
fig.set_size_inches(5, 5)
fig.savefig(str(sc.settings.figdir) + '/umap_lgd_harmony_leiden',
    dpi=400, bbox_extra_artists=(lgd,), bbox_inches='tight')
merge_data.write('harmony_out.h5ad')
