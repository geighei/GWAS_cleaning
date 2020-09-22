#!/usr/bin/env bash
# step 1 - took 5.25 hours on 16cpu-64ram machine; 14 hours on 8cpu-32ram machine
 # /home/ubuntu/tools/regenie/regenie \
 # 	--step 1 \
 # 	--bed ~/UKB/genomeclean/ukb_non_imputed_best_guess_QC \
	# --covarFile ~/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt \
 # 	--phenoFile ~/biroli/geighei/data/GWAS_sumstats/construction/dpw/dpw_phenos.txt \
 # 	--bsize 1000 \
	# --lowmem ~/biroli/geighei/data/GWAS_sumstats/construction/tmpdir/regenie_tmp_preds \
 # 	--out ~/biroli/geighei/data/GWAS_sumstats/construction/dpw/ukb_dpw_test

# step 2 - chromosome 11 took 2 hours on 16cpu-64ram machine; looping chromosomes expected to take ~48 hours
for i in {1..22}
do
	regenie \
	--step 2 \
	--bgen ~/UKB/imputed/ukb_imp_chr${i}_v3.bgen \
	--covarFile ~/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt \
	--phenoFile ~/biroli/geighei/data/GWAS_sumstats/construction/dpw/dpw_phenos.txt \
	--firth 0.01 --approx \
	--bsize 400 \
	--pred ~/biroli/geighei/data/GWAS_sumstats/construction/dpw/full/ukb_dpw_pred.list \
	--split \
	--out ~/biroli/geighei/data/GWAS_sumstats/construction/dpw/full/ukb_dpw_chr${i}
done
