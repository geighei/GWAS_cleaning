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
  - compile_rsid_map.py compiles all RSIDs from the HRC and 1000 genomes panels and combines into a minimal text file for use mapping from CHR:POS coordinates to SNP
  - RS_Merge_Diagnostics.R is a script which provides information about the percentage of SNPs in a sumstat that have RS numbers, by p-value level
  - mergeToRS.py is a program with command-line options which merges the RSIDs legend file (e.g. that from compile_rsid_map.py) to given coordinates of structure CHR:POS to map to RSIDs
  - batch_23andMe.sh is a script that calls mergeToRS.py for all 23&me sumstats
  
  ## SNP_heritability
  - Purpose: Automate computation of SNP heritability of summary statistics using LDSC command-line tool
  
  ## construction
  - Purpose: Construct our own GWASs in UKB using regenie tool and sanity check result with EasyQC
  - preprocess_ukb.py is a script for extracting the necessary components of the UKB phenotypic data for use in UKB GWAS construction in regenie. This mostly means defining and extracting phenotypes and covariates and writing them to individual files
  - single_pheno_regenie.sh is a command-line bash script for running a GWAS on a single phenotype in the UKB using regenie
  - split_sample_regenie.py is a wrapper/master script of single_pheno_regenie.sh. It performs necessary set-up like creating directories. It also allows for running multiple GWASs on an arbitrary number of random "splits" of the sample (e.g. 10-fold cross validation or two splits for ORIV). It also automates the EasyQC diagnostic process after the GWASs are complete.
  - build_ecf.py is a function which reads in a template ECF script and edits the placeholders based on the function inputs and then writes the result to a new ECF script
  - Run_EasyQC.R is a short script that simply runs EasyQC on an ECF script inputted at the command line
  - ecf_scripts/ holds the EasyQC scripts for analysis of the GWASs we construct, including the template script which is used to automatically generate other ECF scripts
