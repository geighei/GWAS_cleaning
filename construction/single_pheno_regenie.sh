#!/usr/bin/env bash
# User-defined inputs
uzh_gwas="/home/ubuntu/biroli/geighei/data/GWAS_sumstats"
bed="/home/ubuntu/UKB/genomeclean/ukb_non_imputed_best_guess_QC"
bgen_prefix="/home/ubuntu/UKB/imputed/ukb_imp_chr"
bgen_suffix="_v3.bgen"
sample_prefix="/home/ubuntu/UKB/genomeraw/ukb41382_imp_chr"
sample_suffix="_v3_s487330.sample"
pheno="dpw"
temp_out="$uzh_gwas/construction/dpw/split_sample/ukb_dpw_splitA"
output="$uzh_gwas/clean/UKB/dpw/ukb_dpw_splitA.txt"

# other stuff
construction="$uzh_gwas/construction"

# step 1 - took 5.25 hours on 16cpu-64ram machine; 14 hours on 8cpu-32ram machine
# extract step: snp list for chr==1-23. added this after new error appeared related to other chr codes
regenie \
--step 1 \
--bed $bed \
--covarFile $construction/ukb_covars.txt \
--phenoFile $construction/$pheno/${pheno}_pheno.txt \
--bsize 1000 \
--lowmem $construction/tmpdir/regenie_tmp_preds \
--out $temp_out

# step 2 - chromosome 11 took 2 hours on 16cpu-64ram machine; looping chromosomes took ~40 hours
for i in {1..22}
do
	regenie \
	--step 2 \
	--bgen ${bgen_prefix}${i}${bgen_suffix} \
	--sample ${sample_prefix}${i}${sample_suffix} \
	--covarFile $construction/ukb_covars.txt \
	--phenoFile $construction/$pheno/${pheno}_pheno.txt \
	--firth 0.01 --approx \
	--bsize 400 \
	--pred ${temp_out}_pred.list \
	--split \
	--out ${temp_out}_chr${i}
done

# step "3" - finalize regenie output for easyQC
# bind together step 2 results, force overwrite
awk 'FNR==1 && NR!=1{next;}{print}' ${temp_out}*${pheno}.regenie > $output

# add N column to resulting file
# store number of samples for which we have the phenotype; this is an estimate, should be good enough
lines=($(wc -l $construction/$pheno/${pheno}_pheno.txt))
# append column with this value in all places to our summary statistics
sed -i '1s/$/ N/; 2,$s/$/'$lines'/' $output
# remove regenie files after doing this because they take up too much space on the drive
rm $temp_out*.regenie