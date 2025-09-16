pyscenic ctx adj_dcm_hcm_healthy.csv \
/data/lixingyu/scenic/vcm/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /data/lixingyu/scenic/vcm/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /data/lixingyu/scenic/vcm/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /data/lixingyu/scenic/vcm/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname /data/lixingyu/scenic/vcm/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /data/lixingyu/scenic/vcm/DCM_HCM_healthy_data.loom \
--output reg_dcm_hcm_healthy.csv \
--mask_dropouts \
--num_workers 10
