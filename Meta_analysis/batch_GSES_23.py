import pandas as pd
from meta23 import meta23
import os
#from shutil import move

# directories
code_dir = "/home/ubuntu/biroli/geighei/code/GWAS_cleaning/Meta_analysis/"
gwas_dir = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/"
gwas_dir_main = gwas_dir + "clean/rGE/"
gwas_dir_23 = gwas_dir + "clean/23andMe/"

# read map file
map = pd.read_csv(gwas_dir + "supplementary/sumstat_maps/cleaned_23andMe_Map.csv")

# loop over rows
for row in map.itertuples(index=False):
	no_23andMe_fp = gwas_dir_main + row.pheno + ".sumstats"
	metal_output_fp = gwas_dir + "clean/metalOutput/" + row.pheno + "_GSES_meta23"
	# if there is no sumstat for file in directory or N is NA, then skip to next
	if (not os.path.isfile(no_23andMe_fp)) | pd.isna(row.N):
		continue

	# Call meta-analysis script using row data as args
	meta23(["--pheno", row.pheno,
	       "--stat1", no_23andMe_fp,
	       "--stat2", gwas_dir_23 + row.file,
	       "--n", str(row.N),
	       "--out", metal_output_fp])
	# Read metal output file
	metal_output = pd.read_csv(metal_output_fp, delim_whitespace=True)
	# Edit column titles to align with PRSice-friendly column titles from main sumstat
	metal_output = metal_output.rename({"MarkerName":"SNP", "Allele1":"A1", "Allele2":"A2",
	                                   "Weight":"N", "Zscore" : "BETA", "P-value" : "P"}, 
	                                   axis="columns")
	# Rename "main" sumstats with _no23andMe suffix
	os.rename(no_23andMe_fp, gwas_dir_main + row.pheno + "_no23andMe.sumstats")
	# Write edited metal output to "main" sumstat place
	metal_output.to_csv(no_23andMe_fp, index=False, sep='\t', na_rep="NA")