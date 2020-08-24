import pandas as pd
import numpy as np
import subprocess
from random import seed

# Steps to 10-fold cross validation w/ 23andMe (#5):
# 1. Read all FID, IID pairs from UKB and randomly label 10 folds
# 2. For each fold:
# a. run regenie (GWAS) on UKB excluding FID, IID pairs in fold
# b. meta-analyze resulting GWAS with _notUKB sumstats (or with both
# 23andMe and _notUKBnot23andMe sumstats)
# c. use resulting sumstats to generate PGS for fold samples
# 3. Append PGS of 10 folds

# set seed so result of GWAS is fully replicable (assignment of folds is random)
seed(112794)

# constant file paths
construction_data = "~/biroli/geighei/data/GWAS_sumstats/construction/"
regenie_exe = "~/tools/regenie/regenie"

# input file paths
covars = construction_data + "ukb_covars.txt"
phenos = construction_data + "dpw_phenos.txt"
geno = "~/UKB/genomeclean/ukb_non_imputed_best_guess_QC"

# read fam file to get list of FID, IID pairs
fam = pd.read_csv(geno + ".fam", 
					names=["FID", "IID", "within_father", "within_mother", "sex", "pheno"], 
					usecols=["FID", "IID"], delim_whitespace=True)
# randomly re-label each sample with unique number and use these to form 10 uniform-sized folds
n_samples = len(fam.index)
fam["random"] = np.random.choice(n_samples, n_samples, replace=False)
fam["fold"] = fam.random % 10

for fold in np.arange(10):
	# select fold subset and write to drive
	fold_ids = fam[fam.fold == fold]
	fold_ids_fp = construction_data + "fold" + str(fold) + ".txt"
	fold_ids.to_csv(fold_ids_fp, sep="\t", index=False)

	# construct regenie calls, excluding the IDs corresponding to this iteration's fold
	step1_cmd = (regenie_exe + " --step 1" + 
				" --remove " + fold_ids_fp +
				" --bed " + geno + 
				" --c " + covars + 
				" --p " + phenos + 
				" --b 1000" + 
				" --o " + construction_data + "dpw_leaveout" + str(fold))
	step2_cmd = (regenie_exe + " --step 2" +
				" --remove " + fold_ids_fp +
				" --bed " + geno + 
				" --c " + covars + 
				" --p " + phenos + 
				" --b 200" +
				" --pred " + construction_data + "dpw_leaveout" + str(fold) + "_pred.list" +
				" --split" +
				" --o ~/biroli/geighei/data/GWAS_sumstats/clean/UKB/dpw_leaveout" + str(fold))

	# call regenie
	subprocess.call(step1_cmd, stderror=subprocess.STDOUT)
	subprocess.call(step2_cmd, stderror=subprocess.STDOUT)

	