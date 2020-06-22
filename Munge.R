##################################
######
###### Project: rGSES
###### 
###### Title: Munge
###### 
###### Author: Michel Navard (edited by Jeremy Vollen)
######
###### Description: Cleans summary statistic files to standardize formatting and, if desired,
###### QC standards. Given a list of files and associated traits, this code cleans, standardizes,
###### and re-writes summary statistic files to the given output directory, logging all changes 
###### to the same directory.
###### 
##################################

# --------------------------------------  Libraries & constants ----------------------------------------------------
library(tidyverse)
library(data.table)
snp_names <- c("snp", "SNP", "snpid", "rsID", "SNPID", "rsid", "RSID", "RS_NUMBER", "rs_number", 
               "RS_NUMBERS", "rs_numbers", "MarkerName", "Markername", "markername", "MARKERNAME", "ID",
               "RS", "SNP_ID", "variant_id")

# --------------------------------------  Helper Functions ----------------------------------------------------

# Reads raw summary statistics and process them (cleanNames) in chunks to avoid using too much memory at one time.
# Calls itself recursively until it detects there is no data left.
# args:
#   read_func: function with which to read summary statistics file, e.g. read_csv, read_table2
#   iter: informs environment how many layers of recursion have transpired, typically default in initial call
#   col_names: informs environment what original column names were so we only read them once, typically default in initial call
#   chunk: how many rows of the data frame to read each iteration
readGWAS <- function(file_path, read_func, iter = 0, col_names = NULL, chunk = 1e7){
  # base case: read first 10 million rows
  if(iter == 0){
    args <- list(file = file_path, col_names=T, na = c(".", NA, "NA", ""), n_max = chunk)
    df <- do.call(read_func, args)
    col_names <- names(df)
  } else { # otherwise, read next 10 million from point we're at, use original col names (passed as arg)
    args <- list(file = file_path, col_names=F, na = c(".", NA, "NA", ""), n_max = chunk, skip = iter*chunk)
    df <- do.call(read_func, args)
    names(df) <- col_names
  }
  # clean this chunk
  df <- cleanNames(df, iter)
  # final case condition, if there weren't chunk rows left, we've reached the end
  if(nrow(df) < chunk) return(df)
  # recursive step. if there's chunk rows, we read the next chunk and bind to this one
  return(bind_rows(df, readGWAS(file_path, read_func, iter = iter + 1, col_names, chunk)))
}

# Function designed for column name exceptions,
# to do general clean up before munging
cleanNames <- function(df, iter) {
  # specifically for Arcogen so we don't have a irregular column name
  names(df) <- str_replace(names(df), "-", "_")
  
  snp_matches <- 
    str_subset(names(df),
               pattern = str_c("^", snp_names, "(_|\\.|$)", collapse = "|"))
  
  rs_found <- 0
  for(match in snp_matches){
    # for each name resembling a SNP column, check if >80% are RS-numbers
    if(sum(str_detect(df[[match]][1:100], "rs"), na.rm=T) > 80){
      # remove all other matching snp names except this one
      df <- df[ , setdiff(names(df), setdiff(snp_matches, match))]
      # rename column containing RS-numbers "SNP"
      data.table::setnames(df, match, "SNP")
      rs_found <- 1
      if(iter == 0){
        # write to log so we know which is the SNP column. invoke calling environment to access log
        cat(print(paste("Interpreting the", match, "column as the SNP column.")),
            file = eval.parent(quote(log.file), 2), sep = "\n", append = TRUE)
      }
      break
    }
  }  
  
  if(!rs_found){
    ## here is where we would write some conversion from CHR:POS to SNP
  }

  return(df)
}


# --------------------------------------  Munge function ----------------------------------------------------

munge <- function (files, trait.names = NULL, out_dir, effect_type = 'both', 
                   N, info.filter = 0, maf.filter = 0) 
  # Input Parameters:
  #   files - list of summary statistic files to be "munged" or cleaned and QC'ed
  #   trait.names - list of phenotype/trait names associated with each summary statistic file in 'files'
  #   out_dir - file path to write "munged" output files and logging file to
  #   effect_type - string describing expected "effect" type; either 'both', 'beta', or 'or'
  #   N - optional argument giving sample size of summary statistic file
  #   info.filter, maf.filter - optional thresholding arguments to perform QC on data
{
  
  ## ----------- Initialize variables and log ---------------
  
  length <- length(files)
  filenames <- as.vector(files)
  
  # Define lists of associated column titles for each desired parameter
  a1_names <- c("a1", "A1", "allele1", "Allele1", "ALLELE1", "EFFECT_ALLELE", "alt", "INC_ALLELE", 
                "REFERENCE_ALLELE", "EA", "Effect_allele", "Effect_Allele", "REF", 
                # added for geighei GWASs
                "effect_allele.exposure", "Tested_Allele", "Allele_1", "effect_allele")
  a2_names <- c("a2", "A2", "allele2", "Allele2", "ALLELE2", "ALLELE0", "OTHER_ALLELE", "ref", 
                "NON_EFFECT_ALLELE", "DEC_ALLELE", "OA", "NEA", "Other_allele", "Other_Allele", 
                "Non_Effect_allele", "Non_Effect_Allele", "ALT", "other_allele.exposure", 
                "Allele_2", "other_allele")
  beta_names <- c("B", "Beta", "beta", "BETA", "EFFECTS", "Effect_Beta", "effect",
                  "EFFECT", "SIGNED_SUMSTAT", "Effect", "b", "est",
                  "Estimate_Effect")
  or_names <- c("OR", "or", "LOG_ODDS", "LogOR", "EFFECTS", "effect", "EFFECT", "SIGNED_SUMSTAT", "Effect", "est")
  # Added to separately account for the common case that a GWAS includes columns for beta
  # and Z. The code suggest renaming Z in this case but we would prefer not to edit the raw data. Instead,
  # we only check for these if no beta name exists
  z_names <- c("Z", "Zscore")
  info_names <-  c("INFO", "info")
  p_names <- c("P", "p", "PVALUE", "Pval", "pvalue", "P_VALUE", "P_value", "P-value", "p-value", 
               "P.value", "p_value", "PVAL", "pval", "P_VAL", "p_val", "GC_PVALUE", "gc_pvalue", 
               "P_Value", "Pvalue")
  n_names <- c("N", "WEIGHT", "nCompleteSamples", "TotalSampleSize", "TotalN", "Total_N", 
               "n_complete_samples", "N_analyzed", "samplesize")
  n_cas_names <- c("NCASE", "N_CASE", "N_CASES", "N_CAS", "NCAS")
  n_con_names <- c("NCONTROL", "N_CONTROL", "N_CONTROLS", "N_CON", "CONTROLS_N", "NCON")
  maf_names <- c("MAF", "maf", "CEUaf", "Freq1", "EAF", "Freq1.Hapmap", "FreqAllele1HapMapCEU", 
                 "Freq.Allele1.HapMapCEU", "EFFECT_ALLELE_FREQ", "Freq.A1")
  
  # Names associated with the effect column are dependent on what effect type the user has requested
  if (effect_type == "beta") {
    effect_names <- beta_names
    effect_name_final <- "BETA"
  } else if (effect_type == "or") {
    effect_names <- or_names
    effect_name_final <- "OR"
  } else if (effect_type == "both") {
    effect_names <- unique(c(beta_names, or_names))
    effect_name_final <- "Effect"
  } else {
    stop("effect_type argument must be either \'both\', \'beta\', or \'or\'")
  }
  
  # Initialize logging file
  log2 <- format(Sys.time(), "%F_%H%M")
  log.file <- file(paste0(out_dir, log2, "_munge.log"), open = "wt")
  begin.time <- Sys.time()
  cat(print(paste0("The munging of ", length(trait.names), 
                   " summary statistics started at ", begin.time), sep = "", quote=F), 
      file = log.file, sep = "\n", append=TRUE)
  
  # Read in each of the files from the inputted list of file paths
  cat(print("Reading and munging the following files:", quote=F), 
      str_c(map(files, print, quote=F)),
      file = log.file, sep = "\n", append = TRUE)

  # Loop through all file/trait names so we generate a summary statistics file for each desired trait
  for (i in 1:length) {
    
    cat(paste("     "), file = log.file, sep = "\n", append = TRUE)
    cat(paste("     "), file = log.file, sep = "\n", append = TRUE)
    cat(print(paste("Munging file:", filenames[i]), quote = F), 
        file = log.file, sep = "\n", append = TRUE)
    
    if (str_sub(files[[i]], start = -3) == "csv") {
      sumstats <- readGWAS(files[[i]], read_csv)
    } else {
      sumstats <- readGWAS(files[[i]], read_table2)
    }
    
    hold_names <- names(sumstats)
    names1 <- hold_names
    
    ## ---------- Find desired columns --------------------------------------  
    # For all columns of interest, check if the desired name already exists;
    #   if so, use that column and write to the log to say you are doing so. 
    # Check if there exist any column names associated with the desired parameter 
    #   and if so, use this column and write to the log what the original name was.
    
    # A1 column
    names1 <- hold_names
    if ("A1" %in% hold_names) 
      cat(print(paste("Interpreting the A1 column as the A1 column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% a1_names] <- "A1"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the", setdiff(names1, 
                                                  hold_names), "column as the A1 column.")), file = log.file, 
          sep = "\n", append = TRUE)
    
    # A2 column
    names1 <- hold_names
    if ("A2" %in% hold_names) 
      cat(print(paste("Interpreting the A2 column as the A2 column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% a2_names] <- "A2"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the", setdiff(names1, 
                                                  hold_names), "column as the A2 column.")), file = log.file, 
          sep = "\n", append = TRUE)
    
    # Effect column
    overlap <- 
      # Flag any column names that are an exact match or start with associated name followed by underscore
      str_detect(hold_names, 
                 pattern = str_c("^", effect_names, "(_|\\.|$)", collapse = "|")) &
      # Remove any column titles containing "allele"; this is necessary because some are called effect_allele_xxxx
      str_detect(hold_names, pattern = "allele", negate = T)
    
    # if the default value already exists, pick it
    if (effect_name_final %in% hold_names) {
      cat(print(paste("Interpreting the", effect_name_final, "column as the effect column.")), 
          file = log.file, sep = "\n", append = TRUE)
      # now ignore case and repeat
    } else if (tolower(effect_name_final) %in% tolower(hold_names)) {
      index <- which(tolower(hold_names) == tolower(effect_name_final))
      cat(print(paste("Interpreting the", hold_names[index], "column as the effect column.")), 
          file = log.file, sep = "\n", append = TRUE)
      hold_names[index] <- effect_name_final
      # else if some other effect name already exists, use that
    } else if (sum(overlap) > 0) {
      cat(print(paste("Interpreting the", 
                      hold_names[overlap], "column as the effect column.")), 
          file = log.file, sep = "\n", append = TRUE)
      hold_names[overlap] <- effect_name_final
      # finally, check if z-score already exists, and if so, use that
    } else if (length(intersect(hold_names, z_names)) > 0) {
      cat(print(paste("Interpreting the", intersect(hold_names,
                                                    z_names), "column as the effect column.")), 
          file = log.file, sep = "\n", append = TRUE)
      hold_names[hold_names %in% z_names] <- effect_name_final
    }
  
    # INFO column
    names1 <- hold_names
    if ("INFO" %in% hold_names) 
      cat(print(paste("Interpreting the INFO column as the INFO column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% info_names] <- "INFO"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the", setdiff(names1, 
                                                  hold_names), "column as the INFO column.")), 
          file = log.file, sep = "\n", append = TRUE)

    # P column - We first check for an exact match, followed by anything starting with p_
    overlap <- str_detect(hold_names, 
                          # Flag any column names that are an exact match or start with associated name followed by underscore
                          pattern = str_c("^", p_names, "(_|\\.|$)", collapse = "|"))
    if ("P" %in% hold_names) 
      cat(print(paste("Interpreting the P column as the P column.")), 
          file = log.file, sep = "\n", append = TRUE)
    else if (sum(overlap) > 0) {
      # sometimes there are numerous reported p-values, in which case we take the linear regression p
      overlap2 <- overlap & str_detect(tolower(hold_names), pattern = "linreg")
      if (sum(overlap) > 1 & sum(overlap2) > 0) {
        overlap <- overlap2
      }
      cat(print(paste("Interpreting the", 
                      hold_names[overlap], "column as the P column.")), 
          file = log.file, sep = "\n", append = TRUE)
      hold_names[overlap] <- "P"
    }
    
    # N column
    names1 <- hold_names
    if ("N" %in% hold_names) 
      cat(print(paste("Interpreting the N column as the N (sample size) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% n_names] <- "N"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the ", setdiff(names1, 
                                                   hold_names), " column as the N (sample size) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # N_CAS column
    names1 <- hold_names
    if ("N_CAS" %in% hold_names) 
      cat(print(paste("Interpreting the N_CAS column as the N_CAS (sample size for cases) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% n_cas_names] <- "N_CAS"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the ", setdiff(names1, 
                                                   hold_names), " column as the N_CAS (sample size for cases) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # N_CON column
    names1 <- hold_names
    if ("N_CON" %in% hold_names) 
      cat(print(paste("Interpreting the N_CON column as the N_CON (sample size for controls) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    hold_names[hold_names %in% n_con_names] <- "N_CON"
    if (length(base::setdiff(names1, hold_names)) > 0) 
      cat(print(paste("Interpreting the ", setdiff(names1, 
                                                   hold_names), " column as the N_CON (sample size for controls) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # If we can't find any column title associated with the desired parameter, then write this to the log
    if (sum(hold_names %in% "P") == 0) 
      cat(print(paste0("Cannot find P-value column, try renaming it to P in the summary statistics file for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "A1") == 0) 
      cat(print(paste0("Cannot find effect allele column, try renaming it to A1 in the summary statistics file for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "A2") == 0) 
      cat(print(paste0("Cannot find other allele column, try renaming it to A2 in the summary statistics file for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% effect_name_final) == 0) 
      cat(print(paste0("Cannot find beta or effect column, try renaming it to effect in the summary statistics file for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "SNP") == 0) 
      cat(print(paste0("Cannot find rs-id column, try renaming it to SNP in the summary statistics file for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    
    # If there wre multiple column titles associated with the desired parameter, suggest column renaming
    if (sum(hold_names %in% "P") > 1) 
      cat(print(paste0("Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "A1") > 1) 
      cat(print(paste0("Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "A2") > 1) 
      cat(print(paste0("Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% effect_name_final) > 1) 
      cat(print(paste0("Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    if (sum(hold_names %in% "SNP") > 1) 
      cat(print(paste0("Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:", 
                       filenames[i])), file = log.file, sep = "\n", 
          append = TRUE)
    
    # Pass warning if we can't find a column associated with one of our desired parameters
    if (sum(hold_names %in% "P") == 0) 
      warning(paste0("Cannot find P-value column, try renaming it P in the summary statistics file for:", 
                     filenames[i]))
    if (sum(hold_names %in% "A1") == 0) 
      warning(paste0("Cannot find effect allele column, try renaming it A1 in the summary statistics file for:", 
                     filenames[i]))
    if (sum(hold_names %in% "A2") == 0) 
      warning(paste0("Cannot find other allele column, try renaming it A2 in the summary statistics file for:", 
                     filenames[i]))
    if (sum(hold_names %in% effect_name_final) == 0) 
      warning(paste0("Cannot find beta or effect column, try renaming it effect in the summary statistics file for:", 
                     filenames[i]))
    if (sum(hold_names %in% "SNP") == 0) 
      warning(paste0("Cannot rs-id column, try renaming it SNP in the summary statistics file for:", 
                     filenames[i]))
    names1 <- hold_names
    
    ## ---------- Quality Control --------------------------------------
    # Edited: similar to effect column, this was made hierarchical because, for example, there are instances
    # where "MAF" and "EAF" occur simultaneously and we of course prefer to use "MAF" and not "EAF"
    if ("MAF" %in% hold_names) { 
      cat(print(paste("Interpreting the MAF column as the MAF (minor allele frequency) column.")), 
          file = log.file, sep = "\n", append = TRUE)
    } else if (length(intersect(hold_names, maf_names)) > 0) {
      cat(print(paste("Interpreting the ", intersect(hold_names, 
                                                     maf_names), " column as the MAF (minor allele frequency) column.")), 
          file = log.file, sep = "\n", append = TRUE)
      hold_names[hold_names %in% maf_names] <- "MAF"
    }

    names(sumstats) <- hold_names
    if ("MAF" %in% colnames(sumstats)) {
      sumstats$MAF <- ifelse(sumstats$MAF <= 0.5, sumstats$MAF, 
                             (1 - sumstats$MAF))
    }
    if ("N_CAS" %in% colnames(sumstats)) {
      sumstats$N <- sumstats$N_CAS + sumstats$N_CON
      cat(print(paste("As the file includes both N_CAS and N_CON columns, the summation of these two columns will be used as the total sample size")), 
          file = log.file, sep = "\n", append = TRUE)
    }
    sumstats$A1 <- factor(toupper(sumstats$A1), c("A", "C", "G", "T"))
    sumstats$A2 <- factor(toupper(sumstats$A2), c("A", "C", "G", "T"))
    
    # Filter out missing p-value rows
    b <- nrow(sumstats)
    if ("P" %in% colnames(sumstats)) {
      sumstats <- subset(sumstats, !(is.na(sumstats$P)))
    }
    if (b - nrow(sumstats) > 0) 
      cat(print(paste(b - nrow(sumstats), "rows were removed from the", 
                      "summary statistics file due to missing values in the P-value column")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # Filter out missing beta value rows
    b <- nrow(sumstats)
    if (effect_name_final %in% colnames(sumstats)) {
      sumstats <- subset(sumstats, !(is.na(sumstats[[effect_name_final]])))
    }
    if (b - nrow(sumstats) > 0) 
      cat(print(paste(b - nrow(sumstats), "rows were removed from the", 
                      "summary statistics file due to missing values in the effect column")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # If beta is determined to be an odds ratio, replace the column with the log of itself
    a1 <- sumstats[[effect_name_final]][[1]]
    sumstats[[effect_name_final]] <- ifelse(rep(round(median(sumstats[[effect_name_final]], 
                                               na.rm = T)) == 1, nrow(sumstats)), log(sumstats[[effect_name_final]]), 
                              sumstats[[effect_name_final]])
    a2 <- sumstats[[effect_name_final]][[1]]
    if (a1 != a2) 
      cat(print(paste("The effect column was determined to be coded as an odds ratio (OR) for the", 
                      filenames[i], "summary statistics file. Please ensure this is correct.
                      The effect column will be of form BETA via computation of log(OR).")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # Sanity check p-value column to make sure no more than 100 values are outside of [0,1] interval
    if ((sum(sumstats$P > 1) + sum(sumstats$P < 0)) > 
        100) {
      cat(print("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!"), 
          file = log.file, sep = "\n", append = TRUE)
    }
    
    # Filter out missing SNP rows
    b <- nrow(sumstats)
    sumstats <- sumstats %>%
      filter(!is.na(SNP) & !is.nan(SNP) & SNP != "NaN")
    if(b - nrow(sumstats) > 0)
      cat(print(paste(b - nrow(sumstats), "rows were removed from the",
                      "summary statistics file due to missing values or NaN in the SNP column")), 
          file = log.file, sep = "\n", append = TRUE)
    
    
    # Filter out any rows with INFO values below user-inputted (or default) threshold
    sumstats$Z <- sign(sumstats[[effect_name_final]]) * sqrt(qchisq(sumstats$P, 1, lower = F))
   
    # Note: we rename OR to BETA at this point since LogOR is treated as a continuous effect
    names(sumstats) <- str_replace_all(names(sumstats), "^OR$", "BETA")
    
    ## ----------- Construct and write output --------------------------------
    # Construct output data frame with relevant columns, including a column for sample size if it doesn't already exist,
    # or pass a warning if we can't find an N column and haven't been passed a sample size
    # if ("N" %in% colnames(sumstats)) {
    output <- sumstats %>%
      as.data.frame() %>%
      select(one_of(c("SNP", "A1", "A2", "BETA", "OR", "Z", "Effect", "P", "CHR", "POS", "N")))

    if (!("N" %in% names(sumstats)) & (exists("N") == FALSE)) 
      cat(warning(paste0("Cannot find sample size column for", 
                         filenames[i], " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.")), 
          file = log.file, sep = "\n", append = TRUE)
    
    # Finalize column names, log about QC changes, and write table to cleaned summary statistics file
    cat(print(paste(nrow(output), "SNPs are left in the summary statistics file", 
                    filenames[i], "after QC.")), file = log.file, sep = "\n", 
        append = TRUE)
    write.table(x = output, file = paste0(out_dir, trait.names[i], ".sumstats"),
                sep = "\t", quote = FALSE, row.names = F)
    #gzip(paste0(trait.names[i], ".sumstats"))
    cat(print(paste("I am done munging file:", filenames[i])), 
        file = log.file, sep = "\n", append = TRUE)
    # cat(print(paste("The file is saved as", paste0(trait.names[i], 
    #                                                ".sumstats.gz"), "in the current working directory.")), 
    #     file = log.file, sep = "\n", append = TRUE)
  }
  
  # Report total time taken and finalize and close the logging file
  end.time <- Sys.time()
  total.time <- difftime(time1 = end.time, time2 = begin.time, 
                         units = "sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time - mins * 60
  cat(paste("     "), file = log.file, sep = "\n", append = TRUE)
  cat(print(paste0("Munging was completed at ", end.time), 
            sep = ""), file = log.file, sep = "\n", append = TRUE)
  cat(print(paste0("The munging of all files took ", mins, 
                   " minutes and ", secs, " seconds"), sep = ""), file = log.file, 
      sep = "\n", append = TRUE)
  cat(print(paste("Please check the log file", paste0(log2, 
                                                      "_munge.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files")), 
      file = log.file, sep = "\n", append = TRUE)
  
  # Close connection to logging file
  close(log.file)
}