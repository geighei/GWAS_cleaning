#!/usr/bin/env bash
shopt -s extglob

##  -------------- Read in options ------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gwas) n="$2"; shift ;;
        --bed) sumstat="$2"; shift ;;
        --bgen1) sumstat="$2"; shift ;;
        --bgen2) sumstat="$2"; shift ;;
        --sample1) sumstat="$2"; shift ;;
        --sample2) sumstat="$2"; shift ;;
        --pheno) sumstat="$2"; shift ;;
        --tmpout) sumstat="$2"; shift ;;
        --out) sumstat="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Example: User-defined inputs for DPW split sample A
# uzh_gwas="/home/ubuntu/biroli/geighei/data/GWAS_sumstats"
# bed="/home/ubuntu/UKB/genomeclean/ukb_non_imputed_best_guess_QC"
# bgen1="/home/ubuntu/UKB/imputed/ukb_imp_chr"
# bgen2="_v3.bgen"
# sample1="/home/ubuntu/UKB/genomeraw/ukb41382_imp_chr"
# sample2="_v3_s487330.sample"
# pheno="dpw"
# tmpout="$uzh_gwas/construction/dpw/split_sample/ukb_dpw_splitA"
# out="$uzh_gwas/clean/UKB/dpw/ukb_dpw_splitA.txt"

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
--out $tmpout

# step 2 - chromosome 11 took 2 hours on 16cpu-64ram machine; looping chromosomes took ~40 hours
for i in {1..22}
do
	regenie \
	--step 2 \
	--bgen ${bgen1}${i}${bgen2} \
	--sample ${sample1}${i}${sample2} \
	--covarFile $construction/ukb_covars.txt \
	--phenoFile $construction/$pheno/${pheno}_pheno.txt \
	--firth 0.01 --approx \
	--bsize 400 \
	--pred ${tmpout}_pred.list \
	--split \
	--out ${tmpout}_chr${i}
done

# step "3" - finalize regenie output for easyQC
# bind together step 2 results, force overwrite
awk 'FNR==1 && NR!=1{next;}{print}' ${tmpout}*${pheno}.regenie > $out

# add N column to resulting file
# store number of samples for which we have the phenotype; this is an estimate, should be good enough
lines=($(wc -l $construction/$pheno/${pheno}_pheno.txt))
# append column with this value in all places to our summary statistics
sed -i '1s/$/ N/; 2,$s/$/'$lines'/' $out
# remove regenie files after doing this because they take up too much space on the drive
rm $tmpout*.regenie