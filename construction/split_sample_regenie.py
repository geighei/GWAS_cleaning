## Name: split_sample_batch.py
##
## Author: Jeremy Vollen
##
## Purpose: automate construction of UKB GWASs using regenie with functionality for 
##		randomly splitting UKB sample into arbitrary number of folds for purpose of 
##		performing Obviously Related Instrumental Variables (ORIV) with resulting PGS
##
## Usage: python split_sample_batch.py -p dpw -n 2 -f split_sample
##		  python split_sample_batch.py -pheno t2d -n 10 -f cv10fold

### IMPORT --------------------------------- ###
import pandas as pd
import numpy as np
import subprocess
import os
import argparse
from contextlib import suppress
from build_ecf import build_ecf

# set seed so result of GWAS is fully replicable (assignment of folds is random)
np.random.seed(112794)

### COMMAND-LINE INPUTS --------------------------------- ###
parser = argparse.ArgumentParser(description='Process inputs.')
parser.add_argument("-p", "--pheno", dest="pheno", required=True,
                    help="name of trait for which we are constructing GWAS")
parser.add_argument("-n", "--n_folds", type=int, dest="n_folds", default=2,
                    help="number of times you want to split UKB for separate GWAS")
parser.add_argument("-f", "--folder", dest="sub_folder", default="split_sample",
                    help="folder name for storing regenie output (goes in phenotype folder in construction)")
args = parser.parse_args()

### FILE PATHS --------------------------------- ###
# constant file paths (note: add 'regenie' executable to your path so you don't need to point to it)
construction_path = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/"
# points to UKB non-imputed data for use in step 1 of regenie
geno_path = "/home/ubuntu/UKB/genomeclean/ukb_non_imputed_best_guess_QC"
# bgen and sample prefix and suffix are the components of the "bgen" and "sample" file paths
# that come before and after the chromosome number
bgen_prefix = "/home/ubuntu/UKB/imputed/ukb_imp_chr"
bgen_suffix = "_v3.bgen"
sample_prefix = "/home/ubuntu/UKB/genomeraw/ukb41382_imp_chr"
sample_suffix = "_v3_s487330.sample"
# filepath where final, binded together regenie output is saved
out_path = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/UKB"

### INTERMEDIATE SETUP --------------------- ###
# input file paths
# phenos file gives phenos but also the selection of individuals to be used in full GWAS
phenos_path = os.path.join(construction_path, args.pheno, args.pheno+"_pheno.txt")
regenie_out_path = os.path.join(construction_path, args.pheno, args.sub_folder)
easyqc_out_path = os.path.join(construction_path, args.pheno, "easyqc")
# make directory for regenie input/output if it doesn't already exist and output directory the same
with suppress(FileExistsError):
	os.makedirs(regenie_out_path)
	os.makedirs(os.path.join(out_path, args.pheno))
	os.makedirs(easyqc_out_path)

# read pheno file to get list of FID, IID pairs that we want to use in analysis
samples = pd.read_csv(phenos_path, delim_whitespace=True)

# randomly re-label each sample with unique number and use these to form n uniform-sized folds
n_samples = len(samples.index)
samples["random"] = np.random.choice(n_samples, n_samples, replace=False)
samples["fold"] = (samples.random % args.n_folds) + 1

### CALL REGENIE --------------------------- ###
for fold in np.arange(1,args.n_folds+1):
	# select fold subset and write to drive
	fold_ids = samples[samples.fold == fold][["FID", "IID"]]
	fold_ids_path = os.path.join(regenie_out_path, "fold"+str(fold)+".txt")
	#fold_ids.to_csv(fold_ids_path, sep="\t", index=False, header=False)
	# save path names which are used more than once
	pheno_split_name = args.pheno+"_split"+str(fold)
	gwas_out_path = os.path.join(out_path, args.pheno, "ukb_"+pheno_split_name+".txt")
	ecf_path = os.path.join("ecf_scripts", "ukb_"+pheno_split_name+".ecf")
	# call script which executes regenie and compiles output according to our specification
	cmd = ["./single_pheno_regenie.sh", 
		"--dir", construction_path,
		"--bed", geno_path,
		"--bgen1", bgen_prefix,
		"--bgen2", bgen_suffix,
		"--sample1", sample_prefix,
		"--sample2", sample_suffix,
		"--pheno", args.pheno,
		"--keep", fold_ids_path,
		"--tmpout", os.path.join(regenie_out_path, "ukb_"+args.pheno),
		"--out", gwas_out_path]
	# run = subprocess.run(cmd)
	
	# RUN EASYQC --------------------------- #
	if fold > 2:
		continue    # For now, we are running EasyQC on first two folds if they exist
	# Build ECF script which will run EasyQC on the GWAS we just ran
	build_ecf(
		template_path=os.path.join("ecf_scripts", "TEMPLATE_easyqc.ecf"),
		out_ecf=ecf_path,
		path_out=easyqc_out_path, 
		file_in=gwas_out_path, 
		sdy=samples.iloc[:,2].std())
	# Call cmd-line R script that calls EasyQC on the ECF script
	run = subprocess.run(["Rscript", "Run_EasyQC.R", ecf_path])