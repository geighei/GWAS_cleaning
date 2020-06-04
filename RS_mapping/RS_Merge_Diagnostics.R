library(readr)
library(tidyverse)
library(stringr)
library(dplyr)
library(purrr)

sumstat_path <- '/home/ubuntu/biroli/geighei/data/GWAS_sumstats/raw/rGE/'

og_filez <- c('osteoarthritis/mapped_GSK_UKBBv3_ArcoGEN-UKHLS.txt',
              'parentalLongevity/deelen/mapped_Results_90th_percentile.txt',
              'type2Diabetes/mahajan/mapped_Mahajan.NatGenet2018b.T2D.European.txt')
pval_names <- c('P-value', 'P-value', 'Pvalue')

evaluateMapping <- function(fp, pval_name, pheno, is_df = F){
  if (is_df)  mapped <- fp 
  else  mapped <- read_table2(str_c(sumstat_path, fp))
  
  names(mapped) <- str_replace_all(names(mapped), pval_name, "pval")
  
  mapped <- mapped %>%
    mutate(unmapped = is.na(SNP),
           top_snp = pval <= 5e-8)
  grouped <- mapped %>%
    group_by(unmapped, top_snp) %>%
    summarise(count = n())
  
  failed_top <- grouped$count[grouped$unmapped == T & grouped$top_snp == T]
  failed_top <- ifelse(length(failed_top) == 0, 0, failed_top)
  failed_other <- grouped$count[grouped$unmapped == T & grouped$top_snp == F]
  
  mapped_top <- grouped$count[grouped$unmapped == F & grouped$top_snp == T]
  mapped_other <- grouped$count[grouped$unmapped == F & grouped$top_snp == F]
  
  total <- nrow(mapped)
  
  writeLines(
    paste(pheno, ":", "\n",
          failed_other + failed_top, "of", total, "total SNPs failed to map to an RS number (",
          (failed_other + failed_top) / total * 100, "% ).\n",
          failed_top, "of", failed_top + mapped_top, "total top SNPs failed to map to an RS number (",
          failed_top / (failed_top + mapped_top) * 100, "% ).\n"))
}

# pl_df <- read_table2(str_c(sumstat_path, og_filez[2]))
# tester <- evaluateMapping(og_filez[2], pval_names[2], 'parentalLongevity')
test <- walk2(og_filez, pval_names,
              ~ evaluateMapping(.x, .y, str_split(.x, "/")[[1]][1]))
