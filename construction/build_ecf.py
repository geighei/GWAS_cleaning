import pandas as pd

def build_ecf(
	# path to template ECF file with placeholders we will replace
	template_path="/home/ubuntu/biroli/geighei/code/GWAS_cleaning/construction/ecf_scripts/TEMPLATE_easyqc.ecf", 
	pheno_path=None, out_ecf, 
	path_out, file_in, sdy=None):
	# error checking
	if pheno_path is None and sdy is None:
		raise ValueError("'build_ecf' function must be supplied either 'pheno_path' or 'sdy' argument")
	# read in file
	with open(template_path, "r") as template:
		template_str = template.read()
	# replace path where EasyQC outputs should be saved
	template_str = template_str.replace("{placeholder_pathOut}", path_out)
	# replace file path to GWAS
	template_str = template_str.replace("{placeholder_fileIn}", file_in)
	# replace file path to GWAS
	template_str = template_str.replace("{placeholder_fileInShortName}", file_in.split("/")[-1].split(".")[0])
	# replace SDY value in string
	if sdy is None:
		pheno = pd.read_csv(pheno_path, delim_whitespace=True)
		sdy = pheno.iloc[:, 2].std()
	template_str = template_str.replace("{placeholder_SDY}", str(sdy))
	# write to file
	with open(out_ecf, 'w') as ecf:
  		ecf.write(template_str)