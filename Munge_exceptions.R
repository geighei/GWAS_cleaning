##################################
######
###### Project: rGSES
###### 
###### Title: Munge Exceptions
###### 
###### Author: Jeremy Vollen
######
###### Description: Some summary statistics defy the logic of munging, usually because of orthogonal standards in 
###### column naming or representing the effect. This file is intended as a place to resolve any GWAS-specific issues that persist
###### after munging. 
###### 
###### 
##################################


flip_list <- c(
  # Liu et al. 2019, alleles must be flipped because beta reported corresponds to alternate rather than reference allele
  'drinkWeek','maxCPD', 'cesSmoke', 'ageSmoke', 'smokeInit', 'dpw')
negate_list <- c(
  # Smoking cessation was defined as 2 for current smokers and 1 for former smokers, the directional opposite of our
  # phenotypic measure, which codes current smokers as 0 and former smokers as 1
  'cesSmoke', 'cancerBreast', 'memoryTest')
master_list <-  c(flip_list, negate_list)

# Read inputted trait's munged sumstat, perform additional operations required for exception, then overwrite sumstat file
postMunge <- function(trait, dir){
  # if trait does not have an exception rule, then exit function without doing anything
  if(trait %not in% master_list){ return() }
  writeLines(str_c("Post: ", trait, " sumstat requires additional cleaning operations"))
  
  sumstat <- read_table2(paste0(dir, trait, '.sumstats'))
  
  if(trait %in% flip_list) {
    sumstat %<>% flipAlleles()
    writeLines(str_c("Post: Flipping alleles of ", trait, " summary statistic"))
  }
  
  if(trait %in% negate_list) {
    sumstat %<>% negateBetas()
    writeLines(str_c("Post: Negating beta of ", trait, " summary statistic"))
  }

  write.table(sumstat, file = paste0(dir, trait, '.sumstats'),
              sep = "\t", quote = FALSE, row.names = F)
  writeLines(str_c("Post: ", trait, " done"))
}

'%not in%' <- function(x,y)!('%in%'(x,y))

# Flip A1 and A2 column, possibly because they were originally misinterpreted in munging process or effect corresponds to A2
flipAlleles <- function(sumstat){
  allele_names <- c('A1', 'A2')
  sumstat[allele_names] <- sumstat[rev(allele_names)]
  return(sumstat)
}

# Negate the effect column, possibly because our phenotype is measured in opposite direction as GWAS sumstat
negateBetas <- function(sumstat){
  sumstat["BETA"] <- -1 * sumstat["BETA"]
  return(sumstat)
}