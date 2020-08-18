#!/usr/bin/env bash

# /home/ubuntu/tools/regenie/regenie \
# 	--step 1 \
# 	--bed ~/UKB/genomeclean/ukb_non_imputed_best_guess_QC \
# 	--c ~/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt \
# 	--p ~/biroli/geighei/data/GWAS_sumstats/construction/dpw_phenos.txt \
# 	--b 1000 \
# 	--o ukb_dpw_test

/home/ubuntu/tools/regenie/regenie \
	--step 2 \
	--bed ~/UKB/genomeclean/ukb_non_imputed_best_guess_QC \
	--c ~/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt \
	--p ~/biroli/geighei/data/GWAS_sumstats/construction/dpw_phenos.txt \
	--b 1000 \
	--pred ukb_dpw_test_pred.list \
	--split \
	--o ukb_dpw_test2
	
