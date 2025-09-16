source activate lxy_py37
cd /Lustre02/lixingyu/test_data/jm_lv
~/anaconda3/envs/lxy_py37/bin/cellbender remove-background \
                 --cuda \
                 --learning-rate 2.5e-5 \
                 --input  /Lustre02/lixingyu/test_data/jm_lv/raw_feature_bc_matrix.h5 \
                 --output /Lustre02/lixingyu/test_data/jm_lv/output.h5
