##################################
######
###### Project: rGSES
###### 
###### Title: GWAS Summary Statistics Cleaning
###### 
###### Author: Jeremy Vollen
######
###### Description: Takes mapping of GWASs to traits and calls munge to produce cleaned
###### and standardized versions of the summary statistics
###### 
##################################


## ----------------------------- Set Path, Macros, Import Libraries -------------------------------
# Clear workspace and load necessary libraries
rm(list = ls())
gc()
library(tidyverse)
library(stringr)
library(magrittr)
library(readxl)
library(optparse)

# Read Command-line arguments
option_list = list(
  # argument pointing to directory where summary statistics live
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="General summary statistics directory", metavar="character"),
  # command-line argument for GWAS file path map
  make_option(c("-m", "--map"), type="character", default=NULL, 
              help="Phenotype filepath map", metavar="character"),
  # argument for directory where we write cleaned GWASs
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output file path", metavar="character"),
  # optional argument: choose subset of phenotypes to clean
  make_option(c("-p", "--pheno_names"), type="character", default=".*",
              help="OPTIONAL: Subset of phenotypes to munge, separated by commas (no spaces), e.g. neuroticismScore,arthritis", 
              metavar="character")) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$map) | is.null(opt$out)) {
  print_help(opt_parser)
  stop("Both a phenotype filepath map and an output file path must be supplied.", call.=FALSE)
}

if (is.null(opt$dir)) {
  print_help(opt_parser)
  stop("Summary statistics file path not provided.")
}

# Load Munge function from separate module
CODE <- "./"
source(paste0(CODE, "Munge.R"))
source(paste0(CODE, "Munge_exceptions.R"))

## ----------------------------- Munge GWASs -------------------------------

# the map provides a mapping from trait (or phenotype) name to summary statistic file path
GWAS_list <- read_csv(opt$map)

# Rename column titles and filter out missing file paths
GWAS_list_clean <- GWAS_list %>%
  # Note: Future versions should force map to give absolute file paths and remove "dir" argument
  select(Pheno_Name, File_Path, Type) %>%
  # Filter out phenotypes for which we don't have sumstat file path
  filter(!is.na(File_Path)) %>%
  # If provided, keep only the user-provided phenotypes; if not, keep all ("." match)
  filter(str_detect(Pheno_Name,
                    str_c("^",
                          str_replace_all(opt$pheno_names, pattern = ",|;", replacement = "$|^"),
                          "$"))) %>%
  # Convert from relative to absolute file path
  mutate(File_Path = str_c(opt$dir, File_Path))

# run OR and BETA GWAS through munge separately since it helps the program discern which effect column to use
beta_GWAS <- GWAS_list_clean %>%
  filter(Type == "beta")
or_GWAS <- GWAS_list_clean %>%
  filter(Type == "or")

# Run munge for all traits requested by the user, running continuous and OR variables in separate calls
if (nrow(beta_GWAS) > 0){
  munge(files = beta_GWAS$File_Path, trait_names = beta_GWAS$Pheno_Name, 
        out_dir = opt$out, effect_type = "beta")
} else { 
  writeLines("Note: None of provided phenotypes are continuous.") 
}

if (nrow(or_GWAS) > 0){
  munge(files = or_GWAS$File_Path, trait_names = or_GWAS$Pheno_Name, 
        out_dir = opt$out, effect_type = "or") 
} else { 
  writeLines("Note: None of provided phenotypes are odds ratios.") 
}

## ----------------------------- Munge exceptions -------------------------------
# run list of traits through munge exception function to clean exception files and overwrite (as side effect, hence 'walk')
GWAS_list_clean$Pheno_Name %>%
  walk(~ postMunge(., dir = opt$out))