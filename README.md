# GWAS_cleaning
All code involving the construction, cleaning, processing, standardizing, "munging", and meta-analysis of GWAS summary statistics.

  ## munging
  - Purpose: standardize/clean sumstats in bulk and write them to files with standardized naming conventions and locations
  - Munge.R stores the 'munge' function, which reads a summary statistic, filters problematic observations, guesses and renames column titles, and writes file to output
  - Munge_exceptions.R performs any common operations that are still required after munging. It isn't very dynamic. It has a list of phenotypes that need to have alleles or beta flipped and does so and overwrites cleaned sumstat files if called on those phenotypes.
  - Clean_GWAS_Cmd.R is a program with command-line options which takes a spreadsheet with sumstat paths and pheno names and munges all phenotypes, writes the standardized versions to a new directory, and calles Munge_exceptions.R for any final touches

  ## meta_analysis
  - Purpose: codes to help automate meta-analysis of two or more GWAS summary statistics
  - meta23.py is a program with command-line options which automates the meta-analysis (using METAL) of a summary statistic with its 23&me cohort sumstat
  - batch_GSES_23.py is a script that runs meta23.py for all phenotypes for which we have 23&me summary statistics
  - metal_commands is a folder where we write text files of commands which are then interpreted by METAl for meta-analysis

  ## RS_mapping
  - Purpose: codes to help mapping RS numbers into summary statistics that don't already have them (e.g. only have CHR-POS coordinates)
  - RS_Merge_Diagnostics.R is a script which provides information about the percentage of SNPs in a sumstat that have RS numbers, by p-value level
  - mergeToRS.py is a program with command-line options which maps RSIDs given coordinates of structure CHR:POS by merging with a combined legend file from the 1000 Genomes project
  - batch_23andMe.sh is a script that calls mergeToRS.py for all 23&me sumstats
  
  ## SNP_heritability
  - Purpose: Automate computation of SNP heritability of summary statistics using LDSC command-line tool
  
  ## construction
  - Purpose: Construct our own GWASs in UKB using regenie tool and sanity check result with EasyQC
  - preprocess_ukb.py is a script for extracting the necessary components of the UKB phenotypic data for use in UKB GWAS construction in regenie. This mostly means defining and extracting phenotypes and covariates and writing them to individual files
  - dpw_ukb_regenie.sh is a script for running a GWAS on drinks per week in the UKB using regenie
  - ten_fold_gwas.py is a script for running a GWAS on a phenotype in UKB using a 10-fold cross-validation technique, necessary if we are to then use this GWAS to construct PGS for individuals in the UKB
  - rGSES_top_phenotypes.xlsx is the spreadsheet where we are currently documenting the UKB variables that map to each selected phenotype for rGSES/GxSES work
  - qc_checks/ holds the EasyQC scripts for analysis of the GWASs we construct
