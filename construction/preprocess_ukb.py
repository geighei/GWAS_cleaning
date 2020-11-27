import pandas as pd
import numpy as np
import re
import os

### DEFINE FILE PATHS
construction_fp = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/construction"
ukb_fp = "/home/ubuntu/biroli/ukb/ukb23283.csv.gz"
sibs_fp = "/home/ubuntu/biroli/ukb/family_siblings_UKB.csv"

#### READ UKB DATA
ukb_cols = ["31",		# Gender
			"34", 		# Year of birth
			"22006",	# Genetic ethnic grouping
			"22027",	# Outliers for heterozygosity or missing rate
			"22019",	# Sex chromosome aneuploidy
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
			"41202", "41204",
						# Type I Diabetes
						# Type II Diabetes
						# Severe Obesity
						# Coronary artery disease (CAD)
						# Alzheimer's Disease
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
			"6150",		# Stroke
			"2405", 	# Number of children fathered (male) 
			"2453",		# Cancer
			"2040",		# Risk taking behaviour
			"41270",	# Alzheimer's
			"6148",		# Cataract
			"2247",		# Hearing difficulty/problems
			"2734",		# Number of live births (female) 
			"20001",	# Number of live births (female) 	
			"20458",	# Positive Affect
			"20460",	# Life Satisfaction	  	
			"20016", 	# Cognitive Performance
			"894", "914", #"884", "904",		
						# Moderate-to-vigorous physical activity
			"2946", "1807", "1845", "3526",
						# Parental longevity
			"20514", "20510", "20517" , "20519", "20511", "20507", "20508", "20518", "20513"      
						# Depressive symptoms
			]

# construct iterator to read zipped file in chunks to minimize computation and memory usage
ukb_iterator = pd.read_csv(ukb_fp, engine="python", encoding = "ISO-8859-1",
					# keep only columns that regex match with our variables of interest since this is a 15GB file
					usecols=lambda col: re.search("eid|^" + "-|^".join(ukb_cols) + "-", col), chunksize=50000)
chunk_list = []
for chunk in ukb_iterator:
	chunk_list.append(chunk)
ukb = pd.concat(chunk_list)

#### CLEAN UKB DATA
# filter on genetically caucasian individuals and filter out heterozygosity, sex outliers
ukb = ukb[(ukb["22006-0.0"] == 1) & (ukb["22027-0.0"] != 1) & (ukb["22019-0.0"] != 1)]
# remove siblings so they can be used as validation
sibs = pd.read_csv(sibs_fp, usecols=["ID"])
non_sibs = set(ukb.eid).difference(sibs.ID)
ukb = ukb[ukb.eid.isin(non_sibs)]

# relabel individual columns
ukb.insert(0, "IID", ukb.eid)
ukb.insert(0, "FID", ukb.eid)
# covariates are gender, year of birth, and first 20 principal components
covar_cols = [col for col in ukb.columns if re.search("^(FID|IID|31-0\.0|34-0\.0|22009-0\.([1-9]$|1[0-9]|20))", col)]

# set of columns to use for all self-reported diagnoses
diagnosis_cols = [col for col in ukb.columns if re.search("^4120(2|4)-", col)]

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
ukb[educ_cols] = ukb[educ_cols].applymap(lambda x: educ_dict.get(x,0))
# take maximum value reported, send to new column
ukb["educYears"] = ukb[educ_cols].max(axis=1)
# filter columns to only keep fid, iid, and education and rows to remove missing education
educ = ukb.dropna(subset=["educYears"])[["FID", "IID", "educYears"]]

# HOUSEHOLD INCOME
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100294
householdIncome_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
householdIncome_cols = [col for col in ukb.columns if re.search("^738-", col)]
ukb[householdIncome_cols] = ukb[householdIncome_cols].applymap(lambda x: householdIncome_dict.get(x,0))
# went with maximum since an average might be skewed by retirement, lay-offs, etc
ukb["householdIncome"] = ukb[householdIncome_cols].max(axis=1)
householdIncome = ukb.dropna(subset=["householdIncome"])[["FID", "IID", "householdIncome"]]

# HEALTH RATING
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100508
health_dict = {1: 4, 2: 3, 3: 2, 4: 1}
health_cols = [col for col in ukb.columns if re.search("^2178-", col)]
ukb[health_cols] = ukb[health_cols].applymap(lambda x: health_dict.get(x,0))
# average health
ukb["healthRating"] = ukb[health_cols].mean(axis=1)
health = ukb.dropna(subset=["healthRating"])[["FID", "IID", "healthRating"]]

# CIGARETTES PER DAY
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2887
# only need to re-code negative values, others are in cigarette units already
maxcpd_dict = {-10: 0, -1: np.nan}
maxcpd_cols = [col for col in ukb.columns if re.search("^2887-", col)]
ukb[maxcpd_cols] = ukb[maxcpd_cols].applymap(lambda x: maxcpd_dict.get(x, x))
# we want to measure propensity to addiction so maximum is more appropriate
ukb["maxcpd"] = ukb[maxcpd_cols].max(axis=1)
maxcpd = ukb.dropna(subset=["maxcpd"])[["FID", "IID", "maxcpd"]]

# AGE FIRST BIRTH (female)
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2754
# only need to re-code negative values, others are in year units already
ageFirstBirth_dict = {-4: np.nan, -3: np.nan}
ageFirstBirth_cols = [col for col in ukb.columns if re.search("^2754-", col)]
ukb[ageFirstBirth_cols] = ukb[ageFirstBirth_cols].applymap(lambda x: ageFirstBirth_dict.get(x, x))
# use first available observation as there shouldn't be inconsistencies
ukb["ageFirstBirth"] = ukb[ageFirstBirth_cols].bfill(axis=1).iloc[:,0]
ageFirstBirth = ukb.dropna(subset=["ageFirstBirth"])[["FID", "IID", "ageFirstBirth"]]

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
ukb["bmi"] = ukb[bmi_cols].mean(axis=1)
bmi = ukb.dropna(subset=["bmi"])[["FID", "IID", "bmi"]]

# SMOKING CESSATION
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=90
cesSmoke_dict = {2: 0, 1: 1}
cesSmoke_cols = [col for col in ukb.columns if re.search("^20116-", col)]
ukb[cesSmoke_cols] = ukb[cesSmoke_cols].applymap(lambda x: cesSmoke_dict.get(x,0))
# use first available observation to maintain consistency across individuals since it's binary
ukb["cesSmoke"] = ukb[cesSmoke_cols].bfill(axis=1).iloc[:,0]
cesSmoke = ukb.dropna(subset=["cesSmoke"])[["FID", "IID", "cesSmoke"]]

# TYPE II DIABETES
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41204
t2d_dict = {"E110": 1, "E111": 1, "E112": 1, "E113": 1, "E114": 1, "E115": 1, "E116": 1, "E117": 1, "E118": 1, "E119": 1}
ukb_t2d = ukb[diagnosis_cols].applymap(lambda x: t2d_dict.get(x, 0))
# select any observation equal to 1
ukb["t2d"] = ukb_t2d.max(axis=1)
t2d = ukb.dropna(subset=["t2d"])[["FID", "IID", "t2d"]]
del ukb_t2d

# TYPE I DIABETES
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41204
t1d_dict = {"E100": 1, "E101": 1, "E102": 1, "E103": 1, "E104": 1, "E105": 1, "E106": 1, "E107": 1, "E108": 1, "E109": 1}
ukb_t1d = ukb[diagnosis_cols].applymap(lambda x: t1d_dict.get(x, 0))
# use first available observation to maintain consistency across individuals since it's binary
ukb["t1d"] = ukb_t1d.max(axis=1)
t1d = ukb.dropna(subset=["t1d"])[["FID", "IID", "t1d"]]
del ukb_t1d

# PROSPECTIVE MEMORY TEST
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20018
memoryTest_dict = {0: 0, 1: 2, 2: 1}
memoryTest_cols = [col for col in ukb.columns if re.search("^20018-", col)]
ukb[memoryTest_cols] = ukb[memoryTest_cols].applymap(lambda x: memoryTest_dict.get(x,0))
# use average value across individuals 
ukb["memoryTest"] = ukb[memoryTest_cols].mean(axis=1)
memoryTest = ukb.dropna(subset=["memoryTest"])[["FID", "IID", "memoryTest"]]

# HIGH BLOOD PRESSURE
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6150
highBloodPressure_dict = {4: 1, -3: np.nan}
highBloodPressure_cols = [col for col in ukb.columns if re.search("^6150-", col)]
ukb[highBloodPressure_cols] = ukb[highBloodPressure_cols].applymap(lambda x: highBloodPressure_dict.get(x,0))
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
ukb[loneliness_cols] = ukb[loneliness_cols].applymap(lambda x: loneliness_dict.get(x,x))
# use first available observation as there shouldn't be inconsistencies
ukb["loneliness"] = ukb[loneliness_cols].max(axis=1)
loneliness = ukb.dropna(subset=["loneliness"])[["FID", "IID", "loneliness"]]

# SMOKE INITIATION
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116
smokeInit_dict = {-3: np.nan, 2: 1}
smokeInit_cols = [col for col in ukb.columns if re.search("^20116-", col)]
ukb[smokeInit_cols] = ukb[smokeInit_cols].applymap(lambda x: smokeInit_dict.get(x,x))
# use last available observation as there shouldn't be inconsistencies
ukb["smokeInit"] = ukb[smokeInit_cols].ffill(axis=1).iloc[:,0]
smokeInit = ukb.dropna(subset=["smokeInit"])[["FID", "IID", "smokeInit"]]

# UNIPOLAR DEPRESSION
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20126
depress_dict = { 0: 0 ,  1: 0 ,  2: 0 , 3: 1,  4: 1, 5: 1 }
depress_cols = [col for col in ukb.columns if re.search("^20126-", col)]
ukb[depress_cols] = ukb[depress_cols].applymap(lambda x: depress_dict.get(x,0))
# use max observation as there shouldn't be inconsistencies
ukb["depress"] = ukb[depress_cols].max(axis=1)
depress = ukb.dropna(subset=["depress"])[["FID", "IID", "depress"]]

# INSOMNIA SYMPTOMS
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1200
insomniaFrequent_dict = {-3: np.nan}
insomniaFrequent_cols = [col for col in ukb.columns if re.search("^1200-", col)]
ukb[insomniaFrequent_cols] = ukb[insomniaFrequent_cols].applymap(lambda x: insomniaFrequent_dict.get(x,x))
# use max observation as there shouldn't be inconsistencies
ukb["insomniaFrequent"] = ukb[insomniaFrequent_cols].max(axis=1)
insomniaFrequent = ukb.dropna(subset=["insomniaFrequent"])[["FID", "IID", "insomniaFrequent"]]

# ARTHRITIS 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20002
arthritis_dict = {1465: 1}
arthritis_cols = [col for col in ukb.columns if re.search("^20002-", col)]
ukb[arthritis_cols] = ukb[arthritis_cols].applymap(lambda x: arthritis_dict.get(x, 0))
# use fist available observation as there shouldn't be inconsistencies
ukb["arthritis"] = ukb[arthritis_cols].max(axis=1)
arthritis = ukb.dropna(subset=["arthritis"])[["FID", "IID", "arthritis"]]

# NON-CANCER ILNESSES
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=135
nonCancerIllness_cols = [col for col in ukb.columns if re.search("^135-", col)]
# use max observation as there shouldn't be inconsistencies
ukb["nonCancerIllness"] = ukb[nonCancerIllness_cols].max(axis=1)
nonCancerIllness = ukb.dropna(subset=["nonCancerIllness"])[["FID", "IID", "nonCancerIllness"]]

# ANXIETY
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20421
anxiety_dict = {-818: np.nan, -121: np.nan}
anxiety_cols = [col for col in ukb.columns if re.search("^20421-", col)]
ukb[anxiety_cols] = ukb[anxiety_cols].applymap(lambda x: anxiety_dict.get(x,x))
# use fist available observation as there shouldn't be inconsistencies
ukb["anxiety"] = ukb[anxiety_cols].max(axis=1)
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
ukb[worryFeeling_cols] = ukb[worryFeeling_cols].applymap(lambda x: worryFeeling_dict.get(x,x))
# use fist available observation as there shouldn't be inconsistencies
ukb["worryFeeling"] = ukb[worryFeeling_cols].max(axis=1)
worryFeeling = ukb.dropna(subset=["worryFeeling"])[["FID", "IID", "worryFeeling"]]

# BREAST CANCER
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=40006
cancerBreast_dict = {"C500": 1, "C501": 1, "C502": 1, "C503": 1, "C504": 1, "C505": 1, "C506": 1, "C507": 1, "C508": 1, "C509": 1}
cancerBreast_cols = [col for col in ukb.columns if re.search("^40006-", col)]
ukb[cancerBreast_cols] = ukb[cancerBreast_cols].applymap(lambda x: cancerBreast_dict.get(x, 0))
# First available observation
ukb["cancerBreast"] = ukb[cancerBreast_cols].max(axis=1)
cancerBreast = ukb.dropna(subset=["cancerBreast"])[["FID", "IID", "cancerBreast"]]

# CHOLESTEROL
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30690
totChol_cols = [col for col in ukb.columns if re.search("^30690-", col)]
# use mean observation as there shouldn't be inconsistencies
ukb["totChol"] = ukb[totChol_cols].mean(axis=1)
totChol = ukb.dropna(subset=["totChol"])[["FID", "IID", "totChol"]]

# STROKE
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6150
stroke_dict = {-3: np.nan, 3: 1}
stroke_cols = [col for col in ukb.columns if re.search("^6150-", col)]
ukb[stroke_cols] = ukb[stroke_cols].applymap(lambda x: stroke_dict.get(x, 0))
# use max observation
ukb["stroke"] = ukb[stroke_cols].max(axis=1)
stroke = ukb.dropna(subset=["stroke"])[["FID", "IID", "stroke"]]

# NUMBER OF CHILDREN FATHERED (MALE) 
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2405
childrenEverFathered_dict = {-1: np.nan, -3: np.nan}
childrenEverFathered_cols = [col for col in ukb.columns if re.search("^2405-", col)]
ukb[childrenEverFathered_cols] = ukb[childrenEverFathered_cols].applymap(lambda x: childrenEverFathered_dict.get(x,x))
# use max observation as there shouldn't be inconsistencies
ukb["childrenEverFathered"] = ukb[childrenEverFathered_cols].max(axis=1)
childrenEverFathered = ukb.dropna(subset=["childrenEverFathered"])[["FID", "IID", "childrenEverFathered"]]

# SEVERE OBESITY
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41204
obesitySevere_dict = {"E660": 1, "E661": 1, "E662": 1, "E663": 1, "E664": 1, "E665": 1, "E666": 1, "E667": 1, "E668": 1, "E669": 1}
ukb_obesitySevere = ukb[diagnosis_cols].applymap(lambda x: obesitySevere_dict.get(x, 0))
# max available observation
ukb["obesitySevere"] = ukb_obesitySevere.max(axis=1)
obesitySevere = ukb.dropna(subset=["obesitySevere"])[["FID", "IID", "obesitySevere"]]
del ukb_obesitySevere

# CANCER
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2453
cancer_dict = {-1: np.nan, -3: np.nan}
cancer_cols = [col for col in ukb.columns if re.search("^2453-", col)]
ukb[cancer_cols] = ukb[cancer_cols].applymap(lambda x: cancer_dict.get(x, x))
# use max observation as there shouldn't be inconsistencies
ukb["cancer"] = ukb[cancer_cols].max(axis=1)
cancer = ukb.dropna(subset=["cancer"])[["FID", "IID", "cancer"]]

# RISK TAKING BEHAVIOUR
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2040
risk_dict = {-1: np.nan, -3: np.nan}
risk_cols = [col for col in ukb.columns if re.search("^2040-", col)]
ukb[risk_cols] = ukb[risk_cols].applymap(lambda x: risk_dict.get(x, x))
# use max observation as there shouldn't be inconsistencies
ukb["risk"] = ukb[risk_cols].max(axis=1)
risk = ukb.dropna(subset=["risk"])[["FID", "IID", "risk"]]

# ALZHEIMER'S
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41270
alzheimer_dict = {"G300": 1, "G301": 1,  "G308": 1, "G309": 1}
ukb_alzheimer = ukb[diagnosis_cols].applymap(lambda x: alzheimer_dict.get(x, 0))
# max available observation
ukb["alzheimer"] = ukb_alzheimer.max(axis=1)
alzheimer = ukb.dropna(subset=["alzheimer"])[["FID", "IID", "alzheimer"]]
del ukb_alzheimer

# CATARACT
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6148
cataract_dict = {-1: np.nan, -3: np.nan, -7: np.nan, 4:1}
# take value 4 , otherwise 0
cataract_cols = [col for col in ukb.columns if re.search("^6148-", col)]
ukb[cataract_cols] = ukb[cataract_cols].applymap(lambda x: cataract_dict.get(x, 0))
# use max observation as there shouldn't be inconsistencies
ukb["cataract"] = ukb[cataract_cols].max(axis=1)
cataract = ukb.dropna(subset=["cataract"])[["FID", "IID", "cataract"]]

# HEARING DIFFICULTY
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2247
hearingDifficulty_dict = {-1: np.nan, -3: np.nan, 99: np.nan}
hearingDifficulty_cols = [col for col in ukb.columns if re.search("^2247-", col)]
ukb[hearingDifficulty_cols] = ukb[hearingDifficulty_cols].applymap(lambda x: hearingDifficulty_dict.get(x,0))
# use max observation as there shouldn't be inconsistencies
ukb["hearingDifficulty"] = ukb[hearingDifficulty_cols].max(axis=1)
hearingDifficulty = ukb.dropna(subset=["hearingDifficulty"])[["FID", "IID", "hearingDifficulty"]]

# NUMBER OF LIVE BIRTH (FEMALE) 
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2734
childrenEverMothered_dict = {-3: np.nan}
childrenEverMothered_cols = [col for col in ukb.columns if re.search("^2734-", col)]
ukb[childrenEverMothered_cols] = ukb[childrenEverMothered_cols].applymap(lambda x: childrenEverMothered_dict.get(x,x))
# use max observation as there shouldn't be inconsistencies
ukb["childrenEverMothered"] = ukb[childrenEverMothered_cols].max(axis=1)
childrenEverMothered = ukb.dropna(subset=["childrenEverMothered"])[["FID", "IID", "childrenEverMothered"]]

# PROSTATE CANCER
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20001
cancerProstate_dict = {1002: 1}
cancerProstate_cols = [col for col in ukb.columns if re.search("^20001-", col)]
ukb[cancerProstate_cols] = ukb[cancerProstate_cols].applymap(lambda x: cancerProstate_dict.get(x, 0))
# use max observation as there shouldn't be inconsistencies
ukb["cancerProstate"] = ukb[cancerProstate_cols].max(axis=1)
cancerProstate = ukb.dropna(subset=["cancerProstate"])[["FID", "IID", "cancerProstate"]]

# CORONARY HEART DISEASE
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41204
cad_dict = {"I250": 1, "I251": 1, "I252": 1, "I253": 1, "I254": 1, "I255": 1, "I256": 1, "I257": 1,  "I258": 1, "I259": 1}
ukb_cad = ukb[diagnosis_cols].applymap(lambda x: cad_dict.get(x, 0))
# max available observation
ukb["cad"] = ukb_cad.max(axis=1)
cad = ukb.dropna(subset=["cad"])[["FID", "IID", "cad"]]
del ukb_cad

# COGNITIVE PERFORMANCE (ALSO DONE BY ANDRIES)
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20016 
cogPerformance_cols = [col for col in ukb.columns if re.search("^20016-", col)]
# use mean cognitive performance as measure 
ukb["cogPerformance"] = ukb[cogPerformance_cols].mean(axis=1)
cogPerformance = ukb.dropna(subset=["cogPerformance"])[["FID", "IID", "cogPerformance"]]

# POSITIVE AFFECT
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20458
positiveAffect_dict = {1:6, 2:5, 3:4, 4:3, 5:2, 6:1}
positiveAffect_cols = [col for col in ukb.columns if re.search("^20458-", col)]
ukb[positiveAffect_cols] = ukb[positiveAffect_cols].applymap(lambda x: positiveAffect_dict.get(x))
# use mean to avoid over-weighting time of particular observation
ukb["positiveAffect"] = ukb[positiveAffect_cols].mean(axis=1)
positiveAffect = ukb.dropna(subset=["positiveAffect"])[["FID", "IID", "positiveAffect"]]

# LIFE SATISFACION
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20460
lifeSatisfaction_dict = {-818: np.nan , -121: np.nan}
lifeSatisfaction_cols = [col for col in ukb.columns if re.search("^20460-", col)]
ukb[lifeSatisfaction_cols] = ukb[lifeSatisfaction_cols].applymap(lambda x: lifeSatisfaction_dict.get(x, x))
# use mean to avoid over-weighting time of particular observation
ukb["lifeSatisfaction"] = ukb[lifeSatisfaction_cols].mean(axis=1)
lifeSatisfaction = ukb.dropna(subset=["lifeSatisfaction"])[["FID", "IID", "lifeSatisfaction"]]

# DEPRESSIVE SYMPTOMS 
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20447
# https://onlinelibrary.wiley.com/doi/full/10.1046/j.1525-1497.2001.016009606.x
depressScore_dict = {1:0, 2:1, 3:2, 4:3}
depressScore_cols = [col for col in ukb.columns if re.search("^(20514|20510|20517|20519|20511|20507|20508|20518|20513)-", col)]
# Subset the columns and standardize so no days of depressive symptomize is 0
ukb_depressScore = ukb[depressScore_cols]
ukb_depressScore = ukb_depressScore.applymap(lambda x: depressScore_dict.get(x))
# Sum all scores 
ukb["depressScore"] = ukb_depressScore.sum(axis=1)
depressScore = ukb.dropna(subset=["depressScore"])[["FID", "IID", "depressScore"]]
del ukb_depressScore

# WELL-BEING SPECTRUM
ukb_wbSpectrum = ukb[["lifeSatisfaction","positiveAffect","neuroticismScore","depressScore"]]
ukb_wbSpectrum_sd = ukb_wbSpectrum.apply(lambda x: (x-x.mean())/x.std())
# would like to use this approach but pandas has an eggregious bug that ignores axis argument
#ukb_wbSpectrum_sd.fillna(value=0, axis=1, limit=1, inplace=True)
ukb_wbSpectrum_sd["nanCount"] = ukb_wbSpectrum_sd.isna().sum(axis=1)
ukb_wbSpectrum_sd.fillna(value=0, inplace=True)
# Take average of normalized wbSpectrum scores, negating neuroticism and depression so they're "positive"
ukb_wbSpectrum_sd["wellBeingSpectrum"] = (
	ukb_wbSpectrum_sd.lifeSatisfaction + ukb_wbSpectrum_sd.positiveAffect - 
	ukb_wbSpectrum_sd.neuroticismScore - ukb_wbSpectrum_sd.depressScore) / 4
ukb["wellBeingSpectrum"] = np.where(ukb_wbSpectrum_sd.nanCount >= 2, np.nan, ukb_wbSpectrum_sd.wellBeingSpectrum)
wellBeingSpectrum = ukb.dropna(subset=["wellBeingSpectrum"])[["FID", "IID", "wellBeingSpectrum"]]
del ukb_wbSpectrum, ukb_wbSpectrum_sd

# AGE PARENTS 90TH
ageParents_cols = [col for col in ukb.columns if re.search("^(2946|1807|1845|3526)-", col)]
ukb_ageParents = ukb[ageParents_cols]
# Combine father and mother's death and age columns to store one value for each
ukb_ageParents["fatherDeath"] = ukb_ageParents.filter(regex="1807").max(axis=1)
ukb_ageParents["motherDeath"] = ukb_ageParents.filter(regex="3526").max(axis=1)
ukb_ageParents["fatherAge"] = ukb_ageParents.filter(regex="2946").max(axis=1)
ukb_ageParents["motherAge"] = ukb_ageParents.filter(regex="1845").max(axis=1)
# Calculate 90th percentile of death age for fathers and mothers separately, excluding Nan's
(father90th, mother90th) = (np.nanpercentile(ukb_ageParents.fatherDeath, 90), np.nanpercentile(ukb_ageParents.motherDeath, 90))
# Individuals marked for parental longevity trait if either mother or father died at or survived to an age >= 90th percentile
ukb_ageParents["ageParents90th"] = 1*(
	(ukb_ageParents.fatherDeath >= father90th) | 
	(ukb_ageParents.motherDeath >= mother90th) |
	(ukb_ageParents.fatherAge >= father90th)   |
	(ukb_ageParents.motherAge >= mother90th))
# NOTE: need to decide still whether we want this (Mark individuals for which we are uncertain as NaNs so they aren't included in GWAS)
ukb["ageParents90th"] = \
	np.where((ukb_ageParents.ageParents90th == 0) & 
		(ukb_ageParents.fatherDeath.isnull() | ukb_ageParents.motherDeath.isnull()), 
		np.nan, ukb_ageParents.ageParents90th)
ageParents = ukb.dropna(subset=["ageParents90th"])[["FID", "IID", "ageParents90th"]]
del ukb_ageParents

# MODERATE TO VIGOROUS PHYSICAL ACTIVITY
actModVig_cols = [col for col in ukb.columns if re.search("^(894|914)-", col)]
# -1 and -3 code missing values; NANs code individuals that reported less than 1 day of 10 min. exercise
ukb_actModVig = ukb[actModVig_cols].fillna(0).replace(to_replace=[-1,-3], value=np.nan)
# For each moderate and vigorous column sets
for reg in ["894", "914"]:
	col_name = "new"+reg
	df = ukb_actModVig.filter(regex=reg)
	# First calculate first non-missing, non-zero value of each row
	col = df.replace(0, np.nan).bfill(1).iloc[:,0]
	# Then, if still missing, set to zero if first value is missing (not -1 or -3)
	# and multiply by 7 to get in min/week instead of min/day
	ukb_actModVig[col_name] = np.where(pd.isnull(col) & (df.iloc[:,0] == 0), 0, 7 * col)
# Weight moderate and vigorous exercise according to MVPA GWAS
ukb["actModVig"] = 4*ukb_actModVig.new894 + 8*ukb_actModVig.new914
actModVig = ukb.dropna(subset=["actModVig"])[["FID", "IID", "actModVig"]]
del ukb_actModVig

# write data
ukb[covar_cols].to_csv(os.path.join(construction_fp, "ukb_covars.txt") sep="\t", index=False, na_rep="NA")
dpw.to_csv(os.path.join(construction_fp, "dpw/dpw_pheno.txt"), sep="\t", index=False, na_rep="NA")
educ.to_csv(os.path.join(construction_fp, "educYears/educYears_pheno.txt"), sep="\t", index=False, na_rep="NA")
householdIncome.to_csv(os.path.join(construction_fp, "householdIncome/householdIncome_pheno.txt"), sep="\t", index=False, na_rep="NA")
health.to_csv(os.path.join(construction_fp, "healthRating/healthRating_pheno.txt"), sep="\t", index=False, na_rep="NA")
maxcpd.to_csv(os.path.join(construction_fp, "maxcpd/maxcpd_pheno.txt"), sep="\t", index=False, na_rep="NA")
ageFirstBirth.to_csv(os.path.join(construction_fp, "ageFirstBirth/ageFirstBirth_pheno.txt"), sep="\t", index=False, na_rep="NA")
smokeInit.to_csv(os.path.join(construction_fp, "smokeInit/smokeInit_pheno.txt"), sep="\t", index=False, na_rep="NA")
bmi.to_csv(os.path.join(construction_fp, "bmi/bmi_pheno.txt"), sep="\t", index=False, na_rep="NA")
cesSmoke.to_csv(os.path.join(construction_fp, "cesSmoke/cesSmoke_pheno.txt"), sep="\t", index=False, na_rep="NA")
t2d.to_csv(os.path.join(construction_fp, "t2d/t2d_pheno.txt"), sep="\t", index=False, na_rep="NA")
t1d.to_csv(os.path.join(construction_fp, "t1d/t1d_pheno.txt"), sep="\t", index=False, na_rep="NA")
memoryTest.to_csv(os.path.join(construction_fp, "memoryTest/memoryTest_pheno.txt"), sep="\t", index=False, na_rep="NA")
highBloodPressure.to_csv(os.path.join(construction_fp, "highBloodPressure/highBloodPressure_pheno.txt"), sep="\t", index=False, na_rep="NA")
medsTaken.to_csv(os.path.join(construction_fp, "medsTaken/medsTaken_pheno.txt"), sep="\t", index=False, na_rep="NA")
loneliness.to_csv(os.path.join(construction_fp, "loneliness/loneliness_pheno.txt"), sep="\t", index=False, na_rep="NA")
smokeInit.to_csv(os.path.join(construction_fp, "smokeInit/smokeInit_pheno.txt"), sep="\t", index=False, na_rep="NA")
depress.to_csv(os.path.join(construction_fp, "depress/depress_pheno.txt"), sep="\t", index=False, na_rep="NA")
insomniaFrequent.to_csv(os.path.join(construction_fp, "insomniaFrequent/insomniaFrequent_pheno.txt"), sep="\t", index=False, na_rep="NA")
arthritis.to_csv(os.path.join(construction_fp, "arthritis/arthritis_pheno.txt"), sep="\t", index=False, na_rep="NA")
nonCancerIllness.to_csv(os.path.join(construction_fp, "nonCancerIllness/nonCancerIllness_pheno.txt"), sep="\t", index=False, na_rep="NA")
anxiety.to_csv(os.path.join(construction_fp, "anxiety/anxiety_pheno.txt"), sep="\t", index=False, na_rep="NA")
height.to_csv(os.path.join(construction_fp, "height/height_pheno.txt"), sep="\t", index=False, na_rep="NA")
asthma.to_csv(os.path.join(construction_fp, "asthma/asthma_pheno.txt"), sep="\t", index=False, na_rep="NA")
neuroticismScore.to_csv(os.path.join(construction_fp, "neuroticismScore/neuroticismScore_pheno.txt"), sep="\t", index=False, na_rep="NA")
worryFeeling.to_csv(os.path.join(construction_fp, "worryFeeling/worryFeeling_pheno.txt"), sep="\t", index=False, na_rep="NA")
cancerBreast.to_csv(os.path.join(construction_fp, "cancerBreast/cancerBreast_pheno.txt"), sep="\t", index=False, na_rep="NA")
totChol.to_csv(os.path.join(construction_fp, "totChol/totChol_pheno.txt"), sep="\t", index=False, na_rep="NA")
stroke.to_csv(os.path.join(construction_fp, "stroke/stroke_pheno.txt"), sep="\t", index=False, na_rep="NA")
childrenEverFathered.to_csv(os.path.join(construction_fp, "childrenEverFathered/childrenEverFathered_pheno.txt"), sep="\t", index=False, na_rep="NA")
obesitySevere.to_csv(os.path.join(construction_fp, "obesitySevere/obesitySevere_pheno.txt"), sep="\t", index=False, na_rep="NA")
cancer.to_csv(os.path.join(construction_fp, "cancer/cancer_pheno.txt"), sep="\t", index=False, na_rep="NA")
risk.to_csv(os.path.join(construction_fp, "risk/risk_pheno.txt"), sep="\t", index=False, na_rep="NA")
alzheimer.to_csv(os.path.join(construction_fp, "alzheimer/alzheimer_pheno.txt"), sep="\t", index=False, na_rep="NA")
cataract.to_csv(os.path.join(construction_fp, "cataract/cataract_pheno.txt"), sep="\t", index=False, na_rep="NA")
hearingDifficulty.to_csv(os.path.join(construction_fp, "hearingDifficulty/hearingDifficulty_pheno.txt"), sep="\t", index=False, na_rep="NA")
childrenEverMothered.to_csv(os.path.join(construction_fp, "childrenEverMothered/childrenEverMothered_pheno.txt"), sep="\t", index=False, na_rep="NA")
cancerProstate.to_csv(os.path.join(construction_fp, "cancerProstate/cancerProstate_pheno.txt"), sep="\t", index=False, na_rep="NA")
cad.to_csv(os.path.join(construction_fp, "cad/cad_pheno.txt"), sep="\t", index=False, na_rep="NA")
cogPerformance.to_csv(os.path.join(construction_fp, "cogPerformance/cogPerformance_pheno.txt"), sep="\t", index=False, na_rep="NA")
positiveAffect.to_csv(os.path.join(construction_fp, "positiveAffect/positiveAffect_pheno.txt"), sep="\t", index=False, na_rep="NA")
lifeSatisfaction.to_csv(os.path.join(construction_fp, "lifeSatisfaction/lifeSatisfaction_pheno.txt"), sep="\t", index=False, na_rep="NA")
depressScore.to_csv(os.path.join(construction_fp, "depressScore/depressScore_pheno.txt"), sep="\t", index=False, na_rep="NA")
wellBeingSpectrum.to_csv(os.path.join(construction_fp, "wellBeingSpectrum/wellBeingSpectrum_pheno.txt"), sep="\t", index=False, na_rep="NA")
ageParents.to_csv(os.path.join(construction_fp, "ageParents90th/ageParents90th_pheno.txt"), sep="\t", index=False, na_rep="NA")
actModVig.to_csv(os.path.join(construction_fp, "actModVig/actModVig_pheno.txt"), sep="\t", index=False, na_rep="NA")