##################################
######
###### Project: rGSES
###### 
###### Title: Summarize UKB phenotypes
###### 
###### Author: Jeremy Vollen
######
###### Description: Compute summary statistics for constructed phenotypes in UKB
###### 
##################################

# Load required packages
requiredPackages = c('tidyverse','magrittr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

options(scipen = 999)
options(digits = 4)
setwd('/Volumes/g_econ_department$/econ/biroli/geighei/data/GWAS_sumstats/construction/')
# number calculated by number of lines in ukb_covars.txt file ("wc -l ukb_covars.txt" in terminal)
total_n = 370008

##
summarizePheno <- function(df){
  pheno <- df[[3]]
  pheno_name <- names(df)[3]
  n <- length(pheno)
  pct_missing <- 1 - (n/total_n)
  mean <- mean(pheno)
  sd <- sd(pheno)
  max <- max(pheno)
  min <- min(pheno)
  output <- tibble(pheno_name = pheno_name, n = n, pct_missing = pct_missing,
                   mean = mean, sd = sd, max = max, min = min)
  return(output)
}


files <- list.files(pattern = "_pheno.txt", recursive = TRUE)
files <- str_subset(files, "cv10fold", negate = TRUE)
df_list <- files %>% 
  map(read_table2)
summary <- df_list %>%
  map(summarizePheno) %>%
  reduce(bind_rows)

write_csv(summary, "ukb_phenos_summary.csv")
