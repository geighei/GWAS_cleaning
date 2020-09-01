import pandas as pd
import numpy as np
import re

#### DPW
# dpw is an exception of sorts since we already had cleaned it some extent in alcohol panel
# panel = pd.read_csv("/home/ubuntu/biroli/geighei/code/alcoholGxE/inputs/alcohol_analysis_panel.csv")
# # include first non-missing wave for each sample, throwing out all obs for which we don't have dpw
# non_missing = panel.dropna(subset=["drinks_per_week"]).groupby(["v_eid"], as_index=False).first()

# non_missing["FID"] = non_missing["v_eid"]
# non_missing["IID"] = non_missing["v_eid"]

# dpw = non_missing[["FID", "IID", "drinks_per_week"]]
# covars = non_missing[["FID", "IID", "male", "year_of_birth", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9", "PC_10",
# 						"PC_11", "PC_12", "PC_13", "PC_14", "PC_15", "PC_16", "PC_17", "PC_18", "PC_19", "PC_20"]]

#### TOP 10 rGSES PHENOTYPES
# del panel, non_missing
ukb_cols = ["eid", 		# Individual ID
			"31_0.0",	# Gender
			"34_0.0", 	# Year of birth
			"22009",	# Principal components
			"6138",		# (Educational) Qualifications
			"738",		# Household income
			"2178",		# Health rating
			"2887",		# Cigs per day
			"2754",		# Age first birth (female)
			"2867",		# Age started smoking
			"21001",	# Body Mass Index (BMI)
			"20116",	# Smoking cessation
			"41204"]	# Type II Diabetes

ukb_iterator = pd.read_csv("/home/ubuntu/biroli/ukb/ukb23283.csv.gz", engine="python", encoding = "ISO-8859-1",
					# keep only columns that regex match with our variables of interest since this is a 15GB file
					usecols=lambda col: re.search("|".join(ukb_cols), col), chunksize=50000)
chunk_list = []
for chunk in ukb_iterator:
	chunk_list.append(chunk)
ukb = pd.concat(chunk_list)

ukb["FID"] = ukb["eid"]
ukb["IID"] = ukb["eid"]
covar_cols = [col for col in ukb.columns if re.search("FID|IID|31_0.0|34_0.0|22009_0_([1-9]$|1[0-9]|20)", col)]


# EDUCATIONAL ATTAINMENT
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100305
# I mapped 5 (NVQ) to nan since it's too ambiguous and 6 is a rough approximation
#educ_dict = {1: 16, 2: 12, 3: 10, 4: 10, 5: np.nan, 6: 14, -7: np.nan, -3: np.nan}
educ_dict = {1: 6, 2: 5, 3: 4, 4: 3, 5: 2, 6: 1, -7: np.nan, -3: np.nan}
# convert all "qualification" columns to years of education
educ_cols = [col for col in ukb.columns if re.search("6138", col)]
ukb[educ_cols] = ukb[educ_cols].applymap(lambda x: educ_dict.get(x))
# take maximum value reported, send to new column
ukb["educYears"] = ukb[educ_cols].max(axis=1)
# filter columns to only keep fid, iid, and education and rows to remove missing education
educ = ukb.dropna(subset=["educYears"])[["FID", "IID", "educYears"]]

# HOUSEHOLD INCOME
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100294
hhi_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
hhi_cols = [col for col in ukb.columns if re.search("738", col)]
ukb[hhi_cols] = ukb[hhi_cols].applymap(lambda x: hhi_dict.get(x))
# went with maximum since an average might be skewed by retirement, lay-offs, etc
ukb["hhi"] = ukb[hhi_cols].max(axis=1)
hhi = ukb.dropna(subset=["hhi"])[["FID", "IID", "hhi"]]

# HEALTH RATING
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100508
health_dict = {1: 4, 2: 3, 3: 2, 4: 1}
health_cols = [col for col in ukb.columns if re.search("2178", col)]
ukb[health_cols] = ukb[health_cols].applymap(lambda x: health_dict.get(x))
# average health
ukb["health_rating"] = ukb[health_cols].mean(axis=1)
health = ukb.dropna(subset=["health_rating"])[["FID", "IID", "health_rating"]]

# CIGARETTES PER DAY
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2887
# only need to re-code negative values, others are in cigarette units already
cpd_dict = {-10: 0, -1: np.nan}
cpd_cols = [col for col in ukb.columns if re.search("2887", col)]
ukb[cpd_cols] = ukb[cpd_cols].applymap(lambda x: cpd_dict.get(x, x))
# we want to measure propensity to addiction so maximum is more appropriate
ukb["cpd"] = ukb[cpd_cols].max(axis=1)
cpd = ukb.dropna(subset=["cpd"])[["FID", "IID", "cpd"]]

# AGE FIRST BIRTH (female)
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2754
# only need to re-code negative values, others are in year units already
afb_dict = {-4: np.nan, -3: np.nan}
afb_cols = [col for col in ukb.columns if re.search("2754", col)]
ukb[afb_cols] = ukb[afb_cols].applymap(lambda x: afb_dict.get(x, x))
# use first available observation as there shouldn't be inconsistencies
ukb["afb"] = ukb[afb_cols].bfill(axis=1).iloc[:,0]
afb = ukb.dropna(subset=["afb"])[["FID", "IID", "afb"]]

# AGE STARTED SMOKING
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2867
# only need to re-code negative values, others are in year units already
smokeInit_dict = {-1: np.nan, -3: np.nan}
smokeInit_cols = [col for col in ukb.columns if re.search("2867", col)]
ukb[smokeInit_cols] = ukb[smokeInit_cols].applymap(lambda x: smokeInit_dict.get(x, x))
# use first available observation as there shouldn't be inconsistencies
ukb["smokeInit"] = ukb[smokeInit_cols].bfill(axis=1).iloc[:,0]
smokeInit = ukb.dropna(subset=["smokeInit"])[["FID", "IID", "smokeInit"]]

# BODY MASS INDEX
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001
bmi_cols = [col for col in ukb.columns if re.search("21001", col)]
# use first available observation, for vast majority this is initial assessment
ukb["bmi"] = ukb[bmi_cols].bfill(axis=1).iloc[:, 0]
bmi = ukb.dropna(subset=["bmi"])[["FID", "IID", "bmi"]]

# SMOKING CESSATION
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=90
cesSmoke_dict = {2: 0, 1: 1}
cesSmoke_cols = [col for col in ukb.columns if re.search("20116", col)]
ukb[cesSmoke_cols] = ukb[cesSmoke_cols].applymap(lambda x: cesSmoke_dict.get(x))
# use first available observation to maintain consistency across individuals since it's binary
ukb["cesSmoke"] = ukb[cesSmoke_cols].bfill(axis=1).iloc[:,0]
cesSmoke = ukb.dropna(subset=["cesSmoke"])[["FID", "IID", "cesSmoke"]]

# TYPE II DIABETES
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41204
t2d_dict = {"E111": 1, "E112": 1, "E113": 1, "E114": 1, "E115": 1, "E116": 1, "E117": 1, "E118": 1, "E119": 1}
t2d_cols = [col for col in ukb.columns if re.search("41204", col)]
ukb[t2d_cols] = ukb[t2d_cols].applymap(lambda x: t2d_dict.get(x, 0))
# use first available observation to maintain consistency across individuals since it's binary
ukb["t2d"] = ukb[t2d_cols].max(axis=1)
t2d = ukb.dropna(subset=["t2d"])[["FID", "IID", "t2d"]]

# write data
# covars.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False, na_rep="NA")
ukb[covar_cols].to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False, na_rep="NA")
# dpw.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/dpw_pheno.txt", sep="\t", index=False, na_rep="NA")
educ.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/educYears_pheno.txt", sep="\t", index=False, na_rep="NA")
hhi.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/householdIncome_pheno.txt", sep="\t", index=False, na_rep="NA")
health.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/healthRating_pheno.txt", sep="\t", index=False, na_rep="NA")
cpd.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/maxCPD_pheno.txt", sep="\t", index=False, na_rep="NA")
afb.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ageFirstBirth_pheno.txt", sep="\t", index=False, na_rep="NA")
smokeInit.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/smokeInit_pheno.txt", sep="\t", index=False, na_rep="NA")
bmi.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/bmi_pheno.txt", sep="\t", index=False, na_rep="NA")
cesSmoke.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/cesSmoke_pheno.txt", sep="\t", index=False, na_rep="NA")
t2d.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/t2d_pheno.txt", sep="\t", index=False, na_rep="NA")


