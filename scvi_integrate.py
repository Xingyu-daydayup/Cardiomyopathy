import pandas as pd
import scanpy as sc
import scvi
annda=sc.read_h5ad('/data/lixingyu/all_data_raw_with_cell_type_full.h5ad')
annda=annda[annda.obs.cell_type1=='vcm']
sample_num=annda.obs.loc[:,'sample'].value_counts()
sample_num=sample_num[sample_num>20]
annda=annda[annda.obs.loc[:,'sample'].isin(sample_num.index)]
annda.layers["counts"] = annda.X.copy()


sc.pp.normalize_total(annda)
sc.pp.log1p(annda)
sc.pp.highly_variable_genes(annda)
annda.raw=annda
annda = annda[:, annda.var.highly_variable]
annda=annda.copy()
scvi.model.SCVI.setup_anndata(annda,layer="counts", batch_key="sample")
model = scvi.model.SCVI(annda, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
annda.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(annda, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(annda)
sc.tl.umap(annda)
annda.write('./vcm_processed.h5ad')




