#!/usr/bin/env bash
shopt -s extglob

# SUMSTATS="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/raw/rGE"
# FILES51="$SUMSTATS/EA3/23andMe/education_years_5.1.dat
# $SUMSTATS/risk/23andMe/risk_preferences.dat
# $SUMSTATS/delayDiscounting/23andMe/MCQ_mean_Log10_K_DecisionMaking_Filtered.5.1.dat
# $SUMSTATS/behaviouralDisinhibition/23andMe/Sanchez-Roige_2017_AUDIT-5.1/AUDIT_Log10_total_clean_DecisionMaking_Filtered.dat"

SUMSTATS="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/raw/23andMe"
FILES51="$SUMSTATS/EA3/education_years_5.1.dat
$SUMSTATS/risk/risk_preferences.dat
$SUMSTATS/delayDiscounting/MCQ_mean_Log10_K_DecisionMaking_Filtered.5.1.dat
$SUMSTATS/behaviouralDisinhibition/AUDIT_Log10_total_clean_DecisionMaking_Filtered.dat"
OUT="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/23andMe"

for file in $FILES51
do python mergeToRS.py --stat $file --map 5.1 --coord all.data.id --is23 --out "$OUT/$(basename "$file")"
done

for file in $SUMSTATS/smoking/!(mapped)*.dat $SUMSTATS/alcohol/!(mapped)*.dat
do python mergeToRS.py --stat $file --map 6.1 --coord all.data.id --is23 --out "$OUT/$(basename "$file")"
done

