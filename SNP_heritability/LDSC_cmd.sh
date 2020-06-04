#!/usr/bin/env bash

# Software: LDSC, information on downloading of the package and reference files can be found here: 
#https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# note python needs to be loaded.

# first a munge step of the sumstats:
gwas="/home/ubuntu/biroli/geighei/data/GWAS_sumstats"
sumstats="$gwas/SNP_heritability/sumstats"
supplementary="$gwas/SNP_heritability/supplementary"
output="$gwas/SNP_heritability/output"
ldsc="/home/ubuntu/plink/ldsc"

# alzheimers
python $ldsc/munge_sumstats.py \
--sumstats $gwas/clean/rGE/alzheimers.sumstats \
--out $sumstats/alzheimers_ldsc \
--signed-sumstats BETA,0 \
--N 455258 \
--merge-alleles $supplementary/w_hm3.noMHC.snplist
# provide additional information on allele names etc if necessary.

# python $ldsc/munge_sumstats.py \
# --sumstats $gwas/clean/rGE/alzheimers.sumstats \
# --out $sumstats/alzheimers_ldsc_altN \
# --signed-sumstats BETA,0 \
# --N 250000 \
# --merge-alleles $supplementary/w_hm3.snplist

## h2 calculation, use output munge as input here:
# Frequency regular drinkers
python $ldsc/ldsc.py \
--h2 $sumstats/alzheimers_ldsc_altN.sumstats.gz \
--ref-ld-chr $supplementary/eur_w_ld_chr/ \
--w-ld-chr $supplementary/eur_w_ld_chr/ \
--out $output/alzheimers

# # Note, for binary traits indicate this make sure to use the liability scale by giving population and sample prevalence, example:
python $ldsc/ldsc.py \
--h2 $sumstats/alzheimers_ldsc_altN.sumstats.gz \
--ref-ld-chr $supplementary/eur_w_ld_chr/ \
--w-ld-chr $supplementary/eur_w_ld_chr/ \
--out $output/alzheimers_binary \
--samp-prev 0.304 \
--pop-prev 0.043



