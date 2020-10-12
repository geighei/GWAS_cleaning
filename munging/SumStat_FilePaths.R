##################################
######
###### Project: rGSES
###### 
###### Title: GWAS Summary Statistics Organization
###### 
###### Author: Jeremy Vollen
######
###### Description: Produces file path list of GWAS summary statistics for all possible phenotypes,
###### which can then be used to more easily map traits to GWASs
###### 
##################################

## ----------------------------- Set Path, Macros, Import Libraries -------------------------------
# Clear workspace and load necessary libraries
rm(list = ls())
gc()
library(tidyverse)
library(stringr)
library(readxl)

# This function reads the column headers from the table at the inputted file path, 
# catching error messages and returning NA
getNames <- function(file_path){
  tryCatch(df <- read_table2(file_path, n_max = 0),
           error = function(cond){
             message(cond)
             return(NA)
           })
  print(file_path)
  return(names(df))
}

# Set working directory to be folder housing GWAS summary statistics
CODE <- "/Volumes/g_econ_department$/econ/biroli/geighei/code/GWAS_cleaning/"
MAP <- "/Volumes/g_econ_department$/econ/biroli/geighei/data/GWAS_sumstats/intermediate/rGE_Phenotype-FilePath_Map.csv"
DATA <- "/Volumes/g_econ_department$/econ/biroli/geighei/data/GWAS_sumstats/raw/"

## ----------------------------- Get all raw files w/ column names -------------------------------

# List all paths in working directory with .txt extension, recursing down file tree
all_txt_files <- list.files(path = paste0(DATA, "rGE"), full.names = T, recursive = T)

# Filter README files out of list and files in OLD directory
no_readme <- all_txt_files[!grepl("readme|old|gz$|png$|docx$|pdf$|xlsx$|tar$|zip$|do$|.me$|pptx$|log$", 
                                  all_txt_files, ignore.case = T)]

# Read and record all column names of each file path remaining, list them in 
# data frame where each row is mapped to a different file path
col_lists <- 
  tibble(file_path = no_readme,
         name_list = map(no_readme, getNames)) %>%
  rowwise() %>%
  mutate(fp_list = list(rep(file_path, times = length(name_list))))

# Unlist the file path and column name lists so the file path is repeated for every 
# column of the .txt file it corresponds to
col_names_df <- tibble(
  file_path = unlist(col_lists$fp_list),
  col_name = unlist(col_lists$name_list))

# Write all raw files that might constitute a GWAS; next step is to manually align GWAS with phenotype of interest
# This is done in the following file: /Users/J_Remedy/Desktop/Career Mode/UZH/Biroli/GWAS Wrangling/Phenotype_FilePath_Map.xlsx
write.table(no_readme, paste0(DATA, "intermediate/all_raw_files.txt"), row.names = F, col.names = F)


## ----------------------------- Pull only those from sumstats we want -------------------------------

# map lists all file paths of sumstats for selected phenotypes/GWASs
setwd(str_c(DATA,".."))
selection <- read_csv(MAP)

# Read and record all column names of selected sumstat file paths,
# list them in df where each row is mapped to a file path
selected_col_lists <- 
  tibble(file_path = selection$File_Path,
         name_list = map(file_path, getNames)) %>%
  rowwise() %>%
  mutate(fp_list = list(rep(file_path, times = length(name_list))))
  
# Unlist the file path and column name lists so the file path is repeated for every 
# column of the .txt file it corresponds to
selected_cols <- 
  tibble(file_path = unlist(selected_col_lists$fp_list),
         col_name = unlist(selected_col_lists$name_list)) %>%
  left_join(selection, by = c("file_path" = "File_Path")) %>%
  select(Phenotype, Pheno_Name, file_path, Type, col_name)

write_csv(selected_cols, "intermediate/selected_cols.csv")
