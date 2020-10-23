import pandas as pd
import numpy as np
import subprocess
import os
from random import seed

# set seed so result of GWAS is fully replicable (assignment of folds is random)
seed(112794)

### INPUTS --------------------------------- ###
n_folds = 2 # this defines the number of times you want to split UKB and GWAS it
# constant file paths (note: add 'regenie' to your path so we don't need to point to it)
construction_fp = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/"
# points to UKB non-imputed data for use in step 1 of regenie
geno_fp = "/home/ubuntu/UKB/genomeclean/ukb_non_imputed_best_guess_QC"
# name of phenotype; eventually, we will extend so we can loop over list of phenotypes
pheno = "dpw"
# sub-folder of this phenotypes folder where we will store results temporarily and write out folds
sub_folder = "split_sample"
out_fp = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/UKB/"dpw/ukb_dpw_splitA.txt"

### INTERMEDIATE SETUP --------------------- ###
# input file paths
# covariates file gives covars but also the selection of individuals to be used in full GWAS
covars_fp = construction_fp + "ukb_covars.txt"
phenos_fp = construction_fp + pheno + "/" + pheno + "_pheno.txt"
regenie_out_fp = construction_fp + pheno + "/" + sub_folder
os.makedirs(regenie_out_fp)

# read fam file to get list of FID, IID pairs
covars = pd.read_csv(covars, 
	usecols=["FID", "IID"], delim_whitespace=True)
# randomly re-label each sample with unique number and use these to form n uniform-sized folds
n_samples = len(covars.index)
covars["random"] = np.random.choice(n_samples, n_samples, replace=False)
covars["fold"] = covars.random % n_folds

### CALL REGENIE --------------------------- ###
for fold in np.arange(1,2): # should be np.arange(1,n_folds+1)
	# select fold subset and write to drive
	fold_ids = covars[covar.fold == fold][["FID", "IID"]]
	fold_ids_fp = regenie_out_fp + "/fold" + str(fold) + ".txt"
	fold_ids.to_csv(fold_ids_fp, sep="\t", index=False, header=False)
	# call script which executes regenie and compiles output according to our specification
	cmd = ["./single_pheno_regenie.sh", 
		"--dir", construction_fp,
		"--bed", fold_ids_fp,
		"--bgen1", "/home/ubuntu/UKB/imputed/ukb_imp_chr",
		"--bgen2", "_v3.bgen",
		"--sample1", "/home/ubuntu/UKB/genomeraw/ukb41382_imp_chr",
		"--sample2", "_v3_s487330.sample",
		"--pheno", pheno,
		"--keep", fold_ids_fp,
		"--tmpout", regenie_out_fp,
		"--out", out_fp + pheno + "/ukb_" + pheno + "_split" + str(fold) + ".txt"]
	run = subprocess.run(cmd)
