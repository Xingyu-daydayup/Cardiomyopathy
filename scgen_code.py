import scgen
import scanpy as sc
annda=sc.read_h5ad('fb_raw.h5ad')
sc.pp.filter_cells(annda, min_genes=200)
sc.pp.filter_genes(annda, min_cells=3)
sc.pp.normalize_total(annda, target_sum=1e4)
sc.pp.log1p(annda)
sc.pp.highly_variable_genes(annda)
annda = annda[:, annda.var.highly_variable]
annda_new=annda.copy()
scgen.setup_anndata(annda_new, batch_key="sample", labels_key="cell_states")
model = scgen.SCGEN(annda_new)
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

corrected_adata = model.batch_removal()
sc.pp.neighbors(corrected_adata)
sc.tl.umap(corrected_adata)
sc.pl.umap(corrected_adata, color=['cell_states'], wspace=0.4, frameon=False,save='test.png')
corrected_adata.write('./fb_data_processed.h5ad')
