import pandas as pd
from meta23 import meta23
import os

# directories
code_dir = "/home/ubuntu/biroli/geighei/code/GWAS_cleaning/Meta_analysis/"
gwas_dir = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/"
gwas_dir_main = gwas_dir + "clean/rGE/"
gwas_dir_23 = gwas_dir + "clean/23andMe/"

# read map file
map = pd.read_csv(gwas_dir + "supplementary/sumstat_maps/cleaned_23andMe_Map.csv")

# loop over rows
for row in map.itertuples(index=False):
	no_23andMe = gwas_dir_main + row.pheno + ".sumstats"
	# if there is no sumstat for file in directory or N is NA, then skip to next
	if (not os.path.isfile(no_23andMe)) | pd.isna(row.N):
		continue
	# call meta-analysis script using row data as args
	meta23(["--pheno", row.pheno,
	       "--stat1", no_23andMe,
	       "--stat2", gwas_dir_23 + row.file,
	       "--n", str(row.N),
	       "--out", gwas_dir + "clean/metaOutput/" + row.pheno + "_GSES_meta23"])
