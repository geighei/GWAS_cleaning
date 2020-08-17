import pandas as pd

panel = pd.read_csv("/home/ubuntu/biroli/geighei/code/alcoholGxE/inputs/alcohol_analysis_panel.csv")

first_wave = panel[panel.i == 1]
first_wave["FID"] = first_wave["v_eid"]
first_wave["IID"] = first_wave["v_eid"]

phenos = first_wave[["FID", "IID", "drinks_per_week"]]
covars = first_wave[["FID", "IID", "male", "year_of_birth", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9", "PC_10"]]

phenos.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/dpw_phenos.txt", sep="\t", index=False)
covars.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False)
