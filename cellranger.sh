#!/bin/bash
/Lustre03/data/wangrui/lixingyu/hospital/cellranger-7.2.0/cellranger count --id=run_jm_lv \
   --fastqs=/Lustre03/data/wangrui/lixingyu/hospital/jm_lv \
   --sample=jm_lv \
   --transcriptome=/Lustre03/data/wangrui/lixingyu/hospital/refdata-gex-GRCh38-2020-A \
   --include-introns true
