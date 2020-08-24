import pandas as pd

panel = pd.read_csv("/home/ubuntu/biroli/geighei/code/alcoholGxE/inputs/alcohol_analysis_panel.csv")
# include first non-missing wave for each sample, throwing out all obs for which we don't have dpw
non_missing = panel.dropna(subset=["drinks_per_week"]).groupby(["v_eid"], as_index=False).first()

non_missing["FID"] = non_missing["v_eid"]
non_missing["IID"] = non_missing["v_eid"]

phenos = non_missing[["FID", "IID", "drinks_per_week"]]
covars = non_missing[["FID", "IID", "male", "year_of_birth", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9", "PC_10"]]

phenos.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/dpw_phenos.txt", sep="\t", index=False, na_rep="NA")
covars.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False, na_rep="NA")
