#!/usr/bin/env bash

construction="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction"
output="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/UKB/dpw"

# step 1 - took 5.25 hours on 16cpu-64ram machine; 14 hours on 8cpu-32ram machine
 # /home/ubuntu/tools/regenie/regenie \
 # 	--step 1 \
 # 	--bed ~/UKB/genomeclean/ukb_non_imputed_best_guess_QC \
  # --covarFile $construction/ukb_covars.txt \
 # 	--phenoFile $construction/dpw/dpw_phenos.txt \
 # 	--bsize 1000 \
	# --lowmem $construction/tmpdir/regenie_tmp_preds \
 # 	--out $construction/dpw/ukb_dpw_test

# step 2 - chromosome 11 took 2 hours on 16cpu-64ram machine; looping chromosomes took ~40 hours
# for i in {1..22}
# do
# 	regenie \
# 	--step 2 \
# 	--bgen ~/UKB/imputed/ukb_imp_chr${i}_v3.bgen \
# 	--sample ~/UKB/genomeraw/ukb41382_imp_chr${i}_v3_s487330.sample \
# 	--covarFile $construction/ukb_covars.txt \
# 	--phenoFile $construction/dpw/dpw_phenos.txt \
# 	--firth 0.01 --approx \
# 	--bsize 400 \
# 	--pred $construction/dpw/full/ukb_dpw_pred.list \
# 	--split \
# 	--out $construction/dpw/full/ukb_dpw_chr${i}
# done

# step "3" - finalize regenie output for easyQC
# bind together step 2 results
awk 'FNR==1 && NR!=1{next;}{print}' $construction/dpw/full/*.regenie > $output/ukb_dpw_full.txt

# add N column to resulting file
	# store number of samples for which we have the phenotype; this is an estimate, should be good enough
lines=($(wc -l $construction/dpw/dpw_phenos.txt))
	# append column with this value in all places to our summary statistics
sed -i '1s/$/ N/; 2,$s/$/'$lines'/' $output/ukb_dpw_full.txt
