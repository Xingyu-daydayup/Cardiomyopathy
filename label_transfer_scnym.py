import scanpy as sc
import numpy as np
from scnym.api import scnym_api
#annda=sc.read_h5ad('../old_data.h5ad')
#annda=annda[annda.obs.cell_type=='fibroblast of cardiac tissue']
annda=sc.read_h5ad('./ref.h5ad')

def subsample(adata, fraction, groupby=None, min_n=0, max_n=10000, method='random', index_only=False, random_state=0):
    if method not in ('random', 'top'):
        raise NotImplementedError(f'method={method} unsupported')
    if groupby:
        if groupby not in adata.obs.columns:
            raise KeyError(f'{groupby} is not a valid obs annotation.')
        groups = adata.obs[groupby].unique()
        n_obs_per_group = {}
        sampled_obs_names = []
        for grp in groups:
            k = adata.obs[groupby] == grp
            grp_size = sum(k)
            ds_grp_size = int(min(
                max_n, max(np.ceil(grp_size * fraction), min(min_n, grp_size))))
            if method == 'top':
                idx = np.argsort(-adata.obs.loc[k, 'n_counts'].values)[0:ds_grp_size]
            else:
                np.random.seed(random_state)
                idx = np.random.choice(grp_size, ds_grp_size, replace=False)
            sampled_obs_names.extend(list(adata.obs_names[k][idx]))
    else:
        ds_size = int(adata.n_obs * fraction)
        np.random.seed(random_state)
        idx = np.random.choice(adata.n_obs, ds_size, replace=False)
        sampled_obs_names = adata.obs_names[idx]
    if index_only:
        return sampled_obs_names
    else:
        return adata[adata.obs_names.isin(sampled_obs_names)].copy()
        
ref = subsample(annda, fraction=1, groupby='cell_states', min_n=0, max_n=2000, method='random', index_only=False, random_state=0)
#marker_genes = [
#    *['MYO18B','XPR1','IGF1R','FNIP2','PCDH7','CCSER1'],
#    *['GPC5','MYL4','TPM3'],
#    *['PRELID2','SH3RF2','GRXCR2'],
#    *['XIRP2','ANKRD2','MYH9'],
#    *['NR4A3','ATF3'],
#    *['CSRP3','ACTB','ATF4'],
#    *['PLXDC2','EBF1','MECOM']
#]
ref.var.index=ref.var.feature_name
#sc.pl.dotplot(ref, marker_genes, 
#              groupby='cell_states',
#              gene_symbols='feature_name',
#              dendrogram=False,
#              # use_raw=True,
#              standard_scale="var",
#              color_map="Reds",
#             swap_axes=False,
#             save='ref.png',
#             show=False)
ref.write('./ref.h5ad')

scnym_api(
    adata=ref,
    task='train',
    groupby='cell_states',
    out_path='./',
    config='new_identity_discovery',
)


query=sc.read_h5ad('./vcm_raw.h5ad')
sc.pp.normalize_total(query, target_sum=1e4)
sc.pp.log1p(query)
scnym_api(
    adata=query,
    task='predict',
    key_added='cell_state_scNym',
    config='new_identity_discovery',
    trained_model='./',
)
#sc.pl.dotplot(query, marker_genes, 
#              groupby='cell_state_scNym',
#              dendrogram=False,
#              # use_raw=True,
#              standard_scale="var",
#              color_map="Reds",
#             swap_axes=False,
#             save='query.png',
#             show=False)
query.write('./vcm_processed.h5ad')
