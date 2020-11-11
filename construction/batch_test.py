import os
import subprocess

#sub_dirs = next(os.walk('.'))[1]
#pheno_dirs = sub_dirs.remove('tmpdirs')
pheno_dirs = ['wellBeingSpectrum', 'educYears', 'householdIncome', 'healthRating', 'maxCPD', 'ageFirstBirth', 'smokeInit', 'bmi', 'cesSmoke', 't2d', 't1d', 'dpw', 'memoryTest', 'highBloodPressure', 'medsTaken', 'loneliness', 'depress', 'insomniaFrequent', 'arthritis', 'nonCancerIllness', 'anxiety', 'height', 'asthma', 'neuroticismScore', 'worryFeeling', 'cancerBreast', 'totChol', 'stroke', 'childrenEverFathered', 'obesitySevere', 'cancer', 'risk', 'alzheimer', 'cataract', 'hearingDifficulty', 'childrenEverMothered', 'cancerProstate', 'cad']

for pheno in pheno_dirs:
	subprocess.call("python split_sample_regenie.py -p " + pheno + " -n 2 -f split_sample")
