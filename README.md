# GWAS_cleaning
All code involving the construction, cleaning, processing, standardizing, "munging", and meta-analysis of GWAS summary statistics.

  ## meta_analysis
  - codes to help automate meta-analysis of two or more GWAS summary statistics
  - meta23.py is a program with command-line options which automates the meta-analysis (using METAL) of a summary statistic with its 23&me cohort sumstat
  - batch_GSES_23.py is a script that runs meta23.py for all phenotypes for which we have 23&me summary statistics
  - metal_commands is a folder where we write text files of commands which are then interpreted by METAl for meta-analysis
