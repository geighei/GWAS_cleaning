##
# Author: Jeremy Vollen
#     Command-line script for calling EasyQC in user-inputted (positional argument)
#		ECF script
##

# Load EasyQC library
library(EasyQC)

# Path to EasyQC ECF script
easyqc_path <- commandArgs(trailingOnly = TRUE)
# Call EasyQC
EasyQC(easyqc_path)