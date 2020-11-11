clear all 

if "`c(username)'" == "debian" {
	global DIR "/home/debian/biroli/geighei/data/GWAS_sumstats/construction"
	}
if "`c(username)'" == "parceh" {
	global DIR  "T:/econ/biroli/geighei/data/GWAS_sumstats/construction"
}


global PHENOS height educYears bmi neuroticismScore worryFeeling risk loneliness insomniaFrequent  depress anxiety  dpw  arthritis t2d ageFirstBirth totChol maxCPD smokeInit cesSmoke cancerBreast cancerProstate  asthma householdIncome cancer cataract cad hearingDifficulty highBloodPressure  medsTaken  childrenEverFathered childrenEverMothered healthRating stroke wellBeingSpectrum t1d  
*actModVig cogPerformance lifeSatisfaction depressScore ageParents90th ageSmoke alzheimer perseveranceLack memoryTest nonCancerIllness

*global PHENOS ageFirstBirth

*global DIR  "T:/econ/biroli/geighei/data/GWAS_sumstats/construction"
*global DIR   "T:\econ\biroli\geighei\data\GWAS_sumstats\construction"

*import delimited using "$DIR/ageFirstBirth/ageFirstBirth_pheno.txt"
set trace off 
foreach i in $PHENOS {
clear
import delimited using "$DIR/`i'/`i'_pheno.txt"

	if "`i'" == "ageFirstBirth" {
	local newvar = lower("`i'")
	rename afb `newvar'
	}
	if "`i'" == "maxCPD" {
	local newvar = lower("`i'")
	rename cpd `newvar'
	
	}
	if "`i'" == "householdIncome" {
	local newvar = lower("`i'")
	rename hhi `newvar'
	
	}
	if "`i'" == "healthRating" {
	local newvar = lower("`i'")
	rename health_rating `newvar'
	
	}
	

local newvar = lower("`i'")
capture confirm variable  `newvar'
if !_rc {
                       di in red "`newvar'  exists"
					   sum  `newvar'
					   clear
               } // end if confirm
               else {
                       di in red "`newvar' does not exist"
					   clear
               } // end else confirm 

	       }
*/
