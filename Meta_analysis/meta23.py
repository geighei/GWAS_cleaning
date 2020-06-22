import pandas as pd
import os
import subprocess

# meta23 - Performs meta-analysis using 23&Me and non-23andMe sumstats by
#			writing METAL command to text file and calling METAL
def meta23(pheno, fp_23, N):
	gwas_fp = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/"
	# open file for writing METAL command script
	cmd_fp = "metal_commands/" + pheno + "_metal.txt"
	cmd = open(cmd_fp, "w")
	# write excl23andMe sumstats file description and process
	cmd.write("MARKER SNP" + "\n" +
	          "DEFAULT " + str(N) + "\n" +
	          "ALLELE A1 A2" + "\n" +
	          "PVAL P" + "\n" +
	          "EFFECT BETA" + "\n" + 
	          "PROCESS " + gwas_fp + "clean/rGE/" + pheno + ".sumstats" + "\n" + "\n")
	# write 23andMe sumstats file description and process
	cmd.write("MARKER SNP" + "\n" +
	          "COLUMN LENIENT" + "\n" +
	          "WEIGHT im.num.0" + "\n" +
	          "ALLELE A1 A2" + "\n" + 
	          "PVAL pvalue" + "\n" +
	          "EFFECT effect" + "\n" +
	          "PROCESS " + gwas_fp + fp_23 + "\n" + "\n")
	# specify output file and analyze
	cmd.write("OUTFILE " + pheno + "_meta_analysis_ " + ".tbl" + "\n" + 
	          "ANALYZE" + "\n")
	cmd.close()

	# run METAL program on script we wrote and output to log
	with open("metal_commands/"+pheno+"_metal.log", "w") as f:
		subprocess.call(["metal", cmd_fp], stdout=f, stderr=subprocess.STDOUT)


N_risk = 466571
meta23("risk", "raw/rGE/risk/23andMe/mapped_risk_preferences.dat", N_risk)