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

#### READ UKB DATA
# del panel, non_missing
ukb_cols = ["eid", 		# Individual ID
			"31-0.0",	# Gender
			"34-0.0", 	# Year of birth
			"22006-0.0",# Genetic ethnic grouping
			"22009",	# Principal components
			"1568", "1578", "1588", "1598", "1608", "5364",
						# Drinks per week
			"6138",		# (Educational) Qualifications
			"738",		# Household income
			"2178",		# Health rating
			"2887",		# Cigs per day
			"2754",		# Age first birth (female)
			"2867",		# Age started smoking
			"21001",	# Body Mass Index (BMI)
			"20116",	# Smoking cessation
			"41204",	# Type II Diabetes
			"20018",	# Prospective memory test 
			"6150",		# High blood pressure
			"137",		# Treatments / medications taken 
			"2020",		# Loneliness 
			"20116",	# Smoke Initiation
			"20126",	# Unipolar Depression
			"1200",		# Insomnia symptoms
  			"20002", 	# Arthritis 
  			"135",		# Non-cancer illnesses 
  			"20421",	# Anxiety 
  			"50",		# Height 
 			"22127",	# Asthma 
 			"20127",	# Neuroticism Score 
			"2000",		# Feeling Worry
			"40006",	# Breast Cancer		
			"30690",	# Cholesterol
			"6150"		# Stroke
]

# construct iterator to read zipped file in chunks to minimize computation and memory usage
ukb_iterator = pd.read_csv("/home/ubuntu/biroli/ukb/ukb23283.csv.gz", engine="python", encoding = "ISO-8859-1",
					# keep only columns that regex match with our variables of interest since this is a 15GB file
					usecols=lambda col: re.search("|".join(ukb_cols), col), chunksize=50000)
chunk_list = []
for chunk in ukb_iterator:
	chunk_list.append(chunk)
ukb = pd.concat(chunk_list)

#### CLEAN UKB DATA
# filter on genetically caucasian individuals and relabel individual columns
ukb = ukb[ukb["22006-0.0"] == 1]
ukb.insert(0, "IID", ukb.eid)
ukb.insert(0, "FID", ukb.eid)
# covariates are gender, year of birth, and first 20 principal components
covar_cols = [col for col in ukb.columns if re.search("^(FID|IID|31-0\.0|34-0\.0|22009-0\.([1-9]$|1[0-9]|20))", col)]
ukb[covar_cols].to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False, na_rep="NA")

# DRINKS PER WEEK
# construction taken from biroli/ukb/alcohol/alcohol_panel_construction/reshape_ukb.R
dpw_dict = {-1: np.nan, -3: np.nan}
dpw_cols = [col for col in ukb.columns if re.search("^(1568|1578|1588|1598|1608|5364)-", col)]
# this is "other beverages" so is more often NA than others and we thus replace NAs with 0
dpw_other_cols = [col for col in ukb.columns if re.search("^5364-", col)]
ukb[dpw_other_cols] = ukb[dpw_other_cols].fillna(0)
# we refer to this several times so let's put it in memory
ukb_dpw = ukb[dpw_cols]
ukb_dpw = ukb_dpw.applymap(lambda x: dpw_dict.get(x, x))
# for each wave, calculate the drinks per week by summing all variables in that wave
for i in range(0,3):
	i = str(i)
	tmp = ukb_dpw.filter(regex="-"+i)
	ukb_dpw["dpw"+i] = tmp.sum(axis=1, skipna=False)
# keep first non-missing wave and put back in ukb dataframe
ukb["dpw"] = ukb_dpw.filter(regex="dpw").bfill(axis=1).iloc[:,0]
dpw = ukb.dropna(subset=["dpw"])[["FID", "IID", "dpw"]]
del ukb_dpw

# EDUCATIONAL ATTAINMENT
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100305
# I mapped 5 (NVQ) to nan since it's too ambiguous and 6 is a rough approximation
#educ_dict = {1: 16, 2: 12, 3: 10, 4: 10, 5: np.nan, 6: 14, -7: np.nan, -3: np.nan}
educ_dict = {1: 6, 2: 5, 3: 4, 4: 3, 5: 2, 6: 1, -7: np.nan, -3: np.nan}
# convert all "qualification" columns to years of education
educ_cols = [col for col in ukb.columns if re.search("^6138-", col)]
ukb[educ_cols] = ukb[educ_cols].applymap(lambda x: educ_dict.get(x))
# take maximum value reported, send to new column
ukb["educYears"] = ukb[educ_cols].max(axis=1)
# filter columns to only keep fid, iid, and education and rows to remove missing education
educ = ukb.dropna(subset=["educYears"])[["FID", "IID", "educYears"]]

# HOUSEHOLD INCOME
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100294
hhi_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
hhi_cols = [col for col in ukb.columns if re.search("^738-", col)]
ukb[hhi_cols] = ukb[hhi_cols].applymap(lambda x: hhi_dict.get(x))
# went with maximum since an average might be skewed by retirement, lay-offs, etc
ukb["hhi"] = ukb[hhi_cols].max(axis=1)
hhi = ukb.dropna(subset=["hhi"])[["FID", "IID", "hhi"]]

# HEALTH RATING
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100508
health_dict = {1: 4, 2: 3, 3: 2, 4: 1}
health_cols = [col for col in ukb.columns if re.search("^2178-", col)]
ukb[health_cols] = ukb[health_cols].applymap(lambda x: health_dict.get(x))
# average health
ukb["health_rating"] = ukb[health_cols].mean(axis=1)
health = ukb.dropna(subset=["health_rating"])[["FID", "IID", "health_rating"]]

# CIGARETTES PER DAY
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2887
# only need to re-code negative values, others are in cigarette units already
cpd_dict = {-10: 0, -1: np.nan}
cpd_cols = [col for col in ukb.columns if re.search("^2887-", col)]
ukb[cpd_cols] = ukb[cpd_cols].applymap(lambda x: cpd_dict.get(x, x))
# we want to measure propensity to addiction so maximum is more appropriate
ukb["cpd"] = ukb[cpd_cols].max(axis=1)
cpd = ukb.dropna(subset=["cpd"])[["FID", "IID", "cpd"]]

# AGE FIRST BIRTH (female)
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2754
# only need to re-code negative values, others are in year units already
afb_dict = {-4: np.nan, -3: np.nan}
afb_cols = [col for col in ukb.columns if re.search("^2754-", col)]
ukb[afb_cols] = ukb[afb_cols].applymap(lambda x: afb_dict.get(x, x))
# use first available observation as there shouldn't be inconsistencies
ukb["afb"] = ukb[afb_cols].bfill(axis=1).iloc[:,0]
afb = ukb.dropna(subset=["afb"])[["FID", "IID", "afb"]]

# AGE STARTED SMOKING
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2867
# only need to re-code negative values, others are in year units already
smokeInit_dict = {-1: np.nan, -3: np.nan}
smokeInit_cols = [col for col in ukb.columns if re.search("^2867-", col)]
ukb[smokeInit_cols] = ukb[smokeInit_cols].applymap(lambda x: smokeInit_dict.get(x, x))
# use first available observation as there shouldn't be inconsistencies
ukb["smokeInit"] = ukb[smokeInit_cols].bfill(axis=1).iloc[:,0]
smokeInit = ukb.dropna(subset=["smokeInit"])[["FID", "IID", "smokeInit"]]

# BODY MASS INDEX
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001
bmi_cols = [col for col in ukb.columns if re.search("^21001-", col)]
# use first available observation, for vast majority this is initial assessment
ukb["bmi"] = ukb[bmi_cols].bfill(axis=1).iloc[:, 0]
bmi = ukb.dropna(subset=["bmi"])[["FID", "IID", "bmi"]]

# SMOKING CESSATION
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=90
cesSmoke_dict = {2: 0, 1: 1}
cesSmoke_cols = [col for col in ukb.columns if re.search("^20116-", col)]
ukb[cesSmoke_cols] = ukb[cesSmoke_cols].applymap(lambda x: cesSmoke_dict.get(x))
# use first available observation to maintain consistency across individuals since it's binary
ukb["cesSmoke"] = ukb[cesSmoke_cols].bfill(axis=1).iloc[:,0]
cesSmoke = ukb.dropna(subset=["cesSmoke"])[["FID", "IID", "cesSmoke"]]

# TYPE II DIABETES
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41204
t2d_dict = {"E110": 1, "E111": 1, "E112": 1, "E113": 1, "E114": 1, "E115": 1, "E116": 1, "E117": 1, "E118": 1, "E119": 1}
t2d_cols = [col for col in ukb.columns if re.search("^41204-", col)]
ukb[t2d_cols] = ukb[t2d_cols].applymap(lambda x: t2d_dict.get(x, 0))
# use first available observation to maintain consistency across individuals since it's binary
ukb["t2d"] = ukb[t2d_cols].max(axis=1)
t2d = ukb.dropna(subset=["t2d"])[["FID", "IID", "t2d"]]

# PROSPECTIVE MEMORY TEST
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20018
memoryTest_dict = {0: 0, 1: 2, 2: 1}
memoryTest_cols = [col for col in ukb.columns if re.search("^20018-", col)]
ukb[memoryTest_cols] = ukb[memoryTest_cols].applymap(lambda x: memoryTest_dict.get(x))
# use average value across individuals 
ukb["memoryTest"] = ukb[memoryTest_cols].mean(axis=1)
memoryTest = ukb.dropna(subset=["memoryTest"])[["FID", "IID", "memoryTest"]]

# HIGH BLOOD PRESSURE
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6150
highBloodPressure_dict = { 1: 0, 2: 0, 3: 0, -3: np.nan, -7: 0}
highBloodPressure_cols = [col for col in ukb.columns if re.search("^6150-", col)]
ukb[highBloodPressure_cols] = ukb[highBloodPressure_cols].applymap(lambda x: highBloodPressure_dict.get(x))
# use maximum to maintain consistency across individuals since it's binary
ukb["highBloodPressure"] = ukb[highBloodPressure_cols].max(axis=1)
highBloodPressure = ukb.dropna(subset=["highBloodPressure"])[["FID", "IID", "highBloodPressure"]]

# TREATMENTS / MADICATIONS TAKEN 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=137
medsTaken_cols = [col for col in ukb.columns if re.search("^137-", col)]
# use average value across individuals 
ukb["medsTaken"] = ukb[medsTaken_cols].mean(axis=1)
medsTaken = ukb.dropna(subset=["medsTaken"])[["FID", "IID", "medsTaken"]]

# LONELINESS
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2020 
loneliness_dict = {-1: np.nan, -3: np.nan}
loneliness_cols = [col for col in ukb.columns if re.search("^2020-", col)]
ukb[loneliness_cols] = ukb[loneliness_cols].applymap(lambda x: loneliness_dict.get(x))
# use first available observation as there shouldn't be inconsistencies
ukb["loneliness"] = ukb[loneliness_cols].bfill(axis=1).iloc[:,0]
loneliness = ukb.dropna(subset=["loneliness"])[["FID", "IID", "loneliness"]]

# SMOKE INITIATION
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116
smokeInit_dict = {-3: np.nan, 2: 1}
smokeInit_cols = [col for col in ukb.columns if re.search("^20116-", col)]
ukb[smokeInit_cols] = ukb[smokeInit_cols].applymap(lambda x: smokeInit_dict.get(x))
# use last available observation as there shouldn't be inconsistencies
ukb["smokeInit"] = ukb[smokeInit_cols].ffill(axis=1).iloc[:,0]
smokeInit = ukb.dropna(subset=["smokeInit"])[["FID", "IID", "smokeInit"]]

# UNIPOLAR DEPRESSION
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20126
depress_dict = { 0: 0 ,  1: 0 ,  2: 0 , 3: 1,  4: 1, 5: 1 }
depress_cols = [col for col in ukb.columns if re.search("^20126-", col)]
ukb[depress_cols] = ukb[depress_cols].applymap(lambda x: depress_dict.get(x))
# use max observation as there shouldn't be inconsistencies
ukb["depress"] = ukb[depress_cols].max(axis=1)
depress = ukb.dropna(subset=["depress"])[["FID", "IID", "depress"]]

# INSOMNIA SYMPTOMS
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1200
insomniaFrequent_dict = {-3: np.nan}
insomniaFrequent_cols = [col for col in ukb.columns if re.search("^1200-", col)]
ukb[insomniaFrequent_cols] = ukb[insomniaFrequent_cols].applymap(lambda x: insomniaFrequent_dict.get(x))
# use max observation as there shouldn't be inconsistencies
ukb["insomniaFrequent"] = ukb[insomniaFrequent_cols].max(axis=1)
insomniaFrequent = ukb.dropna(subset=["insomniaFrequent"])[["FID", "IID", "insomniaFrequent"]]

# ARTHRITIS 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20002
arthritis_dict = {99999: np.nan, 1465: 1}
arthritis_cols = [col for col in ukb.columns if re.search("^20002-", col)]
ukb[arthritis_cols] = ukb[arthritis_cols].applymap(lambda x: arthritis_dict.get(x, x))
# use fist available observation as there shouldn't be inconsistencies
ukb["arthritis"] = ukb[arthritis_cols].bfill(axis=1).iloc[:,0]
arthritis = ukb.dropna(subset=["arthritis"])[["FID", "IID", "arthritis"]]

# NON-CANCER ILNESSES
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=135
nonCancerIllness_cols = [col for col in ukb.columns if re.search("^135-", col)]
# use max observation as there shouldn't be inconsistencies
ukb["nonCancerIllness"] = ukb[nonCancerIllness_cols].max(axis=1)
nonCancerIllness_cols = ukb.dropna(subset=["nonCancerIllness"])[["FID", "IID", "nonCancerIllness"]]

# ANXIETY
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20421
anxiety_dict = {-818: np.nan, -121: np.nan}
anxiety_cols = [col for col in ukb.columns if re.search("^20421-", col)]
ukb[anxiety_cols] = ukb[anxiety_cols].applymap(lambda x: anxiety_dict.get(x))
# use fist available observation as there shouldn't be inconsistencies
ukb["anxiety"] = ukb[anxiety_cols].bfill(axis=1).iloc[:,0]
anxiety = ukb.dropna(subset=["anxiety"])[["FID", "IID", "anxiety"]]

# HEIGHT
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=50
height_cols = [col for col in ukb.columns if re.search("^50-", col)]
# use max observation as there shouldn't be inconsistencies
ukb["height"] = ukb[height_cols].max(axis=1)
height = ukb.dropna(subset=["height"])[["FID", "IID", "height"]]

# ASTHMA
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22127
asthma_cols = [col for col in ukb.columns if re.search("^22127-", col)]
# use max observation as there shouldn't be inconsistencies
ukb["asthma"] = ukb[asthma_cols].max(axis=1)
asthma = ukb.dropna(subset=["asthma"])[["FID", "IID", "asthma"]]

# NEUROTICISM SCORE
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20127
neuroticismScore_cols = [col for col in ukb.columns if re.search("^20127-", col)]
# use max observation as there shouldn't be inconsistencies
ukb["neuroticismScore"] = ukb[neuroticismScore_cols].max(axis=1)
neuroticismScore = ukb.dropna(subset=["neuroticismScore"])[["FID", "IID", "neuroticismScore"]]

# FEELING WORRY
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2000
worryFeeling_dict = {-1: np.nan, -3: np.nan}
worryFeeling_cols = [col for col in ukb.columns if re.search("^2000-", col)]
ukb[worryFeeling_cols] = ukb[worryFeeling_cols].applymap(lambda x: worryFeeling_dict.get(x))
# use fist available observation as there shouldn't be inconsistencies
ukb["worryFeeling"] = ukb[worryFeeling_cols].max(axis=1)
worryFeeling = ukb.dropna(subset=["worryFeeling"])[["FID", "IID", "worryFeeling"]]

# BREAST CANCER
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=40006
cancerBreast_dict = {"C500": 1, "C501": 1, "C502": 1, "C503": 1, "C504": 1, "C505": 1, "C506": 1, "C507": 1, "C508": 1, "C509": 1}
cancerBreast_cols = [col for col in ukb.columns if re.search("^40006-", col)]
ukb[cancerBreast_cols] = ukb[cancerBreast_cols].applymap(lambda x: cancerBreast_dict.get(x, 0))
# First available observation
ukb["cancerBreast"] = ukb[cancerBreast_cols].bfill(axis=1).iloc[:,0]
cancerBreast = ukb.dropna(subset=["cancerBreast"])[["FID", "IID", "cancerBreast"]]

# CHOLESTEROL
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30690
totChol_cols = [col for col in ukb.columns if re.search("^30690-", col)]
# use mean observation as there shouldn't be inconsistencies
ukb["totChol"] = ukb[totChol_cols].mean(axis=1)
totChol = ukb.dropna(subset=["totChol"])[["FID", "IID", "totChol"]]

# STROKE
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6150
stroke_dict = {1:0, 3:1}
stroke_cols = [col for col in ukb.columns if re.search("^6150-", col)]
ukb[stroke_cols] = ukb[stroke_cols].applymap(lambda x: stroke_dict.get(x, 0))
# use first available observation
ukb["stroke"] = ukb[stroke_cols].bfill(axis=1).iloc[:,0]
stroke = ukb.dropna(subset=["stroke"])[["FID", "IID", "stroke"]]

# write data
ukb[covar_cols].to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ukb_covars.txt", sep="\t", index=False, na_rep="NA")
dpw.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/dpw_pheno.txt", sep="\t", index=False, na_rep="NA")
educ.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/educYears_pheno.txt", sep="\t", index=False, na_rep="NA")
hhi.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/householdIncome_pheno.txt", sep="\t", index=False, na_rep="NA")
health.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/healthRating_pheno.txt", sep="\t", index=False, na_rep="NA")
cpd.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/maxCPD_pheno.txt", sep="\t", index=False, na_rep="NA")
afb.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/ageFirstBirth_pheno.txt", sep="\t", index=False, na_rep="NA")
smokeInit.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/smokeInit_pheno.txt", sep="\t", index=False, na_rep="NA")
bmi.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/bmi_pheno.txt", sep="\t", index=False, na_rep="NA")
cesSmoke.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/cesSmoke_pheno.txt", sep="\t", index=False, na_rep="NA")
t2d.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/t2d_pheno.txt", sep="\t", index=False, na_rep="NA")
memoryTest.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/memoryTest_pheno.txt", sep="\t", index=False, na_rep="NA")
highBloodPressure.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/highBloodPressure_pheno.txt", sep="\t", index=False, na_rep="NA")
medsTaken.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/medsTaken_pheno.txt", sep="\t", index=False, na_rep="NA")
loneliness.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/loneliness_pheno.txt", sep="\t", index=False, na_rep="NA")
smokeInit.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/smokeInit_pheno.txt", sep="\t", index=False, na_rep="NA")
depress.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/depress_pheno.txt", sep="\t", index=False, na_rep="NA")
insomniaFrequent.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/insomniaFrequent_pheno.txt", sep="\t", index=False, na_rep="NA")
arthritis.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/arthritis_pheno.txt", sep="\t", index=False, na_rep="NA")
nonCancerIllness.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/nonCancerIllness_pheno.txt", sep="\t", index=False, na_rep="NA")
anxiety.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/anxiety_pheno.txt", sep="\t", index=False, na_rep="NA")
height.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/height_pheno.txt", sep="\t", index=False, na_rep="NA")
asthma.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/asthma_pheno.txt", sep="\t", index=False, na_rep="NA")
neuroticismScore.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/neuroticismScore_pheno.txt", sep="\t", index=False, na_rep="NA")
worryFeeling.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/worryFeeling_pheno.txt", sep="\t", index=False, na_rep="NA")
cancerBreast.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/cancerBreast_pheno.txt", sep="\t", index=False, na_rep="NA")
totChol.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/totChol_pheno.txt", sep="\t", index=False, na_rep="NA")
stroke.to_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction/stroke_pheno.txt", sep="\t", index=False, na_rep="NA")
