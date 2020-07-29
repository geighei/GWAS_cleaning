import pandas as pd
import os
import subprocess
import argparse

## INPUTS: file paths to summary statistic, 1000G legend, and desired output file path
def meta23(raw_args=None):
     parser = argparse.ArgumentParser(description='Process file paths and relevant column names for merging.')
     parser.add_argument("-p", "--pheno", dest="pheno", required=True,
                         help="name of trait for which we are performing meta-analysis")
     parser.add_argument("--stat1", dest="fp_no23",
                         help="file path to first summary statistic file, usually original")
     parser.add_argument("--stat2", dest="fp_23",
                         help="file path to second summary statistic file, usually component missing from GWAS")
     parser.add_argument("-n", "--n", dest="N", default=100000,
                         help="number of samples included in summary statistic #1")
     parser.add_argument("-o", "--out", dest="output_fp", default = "./meta_output",
                         help="file path for meta-analysis output")
     args = parser.parse_args(raw_args)

     # meta23 - Performs meta-analysis using 23&Me and non-23andMe sumstats by
     #			writing METAL command to text file and calling METAL
     code_fp = "/home/ubuntu/biroli/geighei/code/GWAS_cleaning/Meta_analysis/"
     gwas_fp = "/home/ubuntu/biroli/geighei/data/GWAS_sumstats/"
     # open file for writing METAL command script
     cmd_fp = code_fp + "metal_commands/" + args.pheno + "_metal.txt"
     cmd = open(cmd_fp, "w")
     # write excl23andMe sumstats file description and process
     cmd.write("MARKER SNP" + "\n" +
               "DEFAULT " + str(args.N) + "\n" +
               "WEIGHT N" + "\n" +
               "ALLELE A1 A2" + "\n" +
               "PVAL P" + "\n" +
               "EFFECT BETA" + "\n" + 
               "PROCESS " + args.fp_no23 + "\n" + "\n")
     # write 23andMe sumstats file description and process
     cmd.write("MARKER SNP" + "\n" +
               "COLUMN LENIENT" + "\n" +
               "WEIGHT im.num.0" + "\n" +
               "ALLELE A1 A2" + "\n" + 
               "PVAL pvalue" + "\n" +
               "EFFECT effect" + "\n" +
               "PROCESS " + args.fp_23 + "\n" + "\n")
     # specify output file and analyze
     cmd.write("OUTFILE " + args.output_fp + "_ .tbl" + "\n" + 
               "ANALYZE" + "\n")
     cmd.close()

     # run METAL program on script we wrote and output to log
     with open(code_fp+"metal_commands/"+args.pheno+"_metal.log", "w") as f:
     	subprocess.call(["metal", cmd_fp], stdout=f, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    meta23()