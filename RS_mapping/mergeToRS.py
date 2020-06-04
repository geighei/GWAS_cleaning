#### This script maps RSIDs given coordinates of structure CHR:POS
####   To do this, it merges the given summary statistic path with 
####   a combined legend file from the 1000 Genomes project

import pandas as pd
import numpy as np 
import os
import argparse

## INPUTS: file paths to summary statistic, 1000G legend, and desired output file path
parser = argparse.ArgumentParser(description='Process file paths and relevant column names for merging.')
parser.add_argument("-s", "--stat", dest="sumstat_fp", required=True,
                    help="file path to summary statistic file to be munged")
parser.add_argument("-m", "--map", dest="legend_fp", default="/home/ubuntu/dep/rsid_map/rsid_map.txt",
                    help="file path to CHR:POS -> RSID map; if 23andMe, just provide version number (e.g. 5.1)")
parser.add_argument("-c", "--coord", dest="sumstat_coord_name", default="SNP",
                    help="name of column in summary statistic file which corresponds to CHR:POS coordinates")
parser.add_argument("--is23", dest="is_23andMe", action='store_true',
                    help="boolean for whether mapping is for 23andMe-provided sumstat")
parser.add_argument("-o", "--out", dest="output_fp", default = None,
                    help="file path for mapped summary statistic file")

args = parser.parse_args()

## ARGUMENT CLEANING AND EXCEPTIONS
# if output file path is not supplied, write to same location as sumstat but with prefix
if args.output_fp == None:
	args.output_fp = (os.path.dirname(args.sumstat_fp) + 
	                  "/mapped_" + os.path.basename(args.sumstat_fp))

# if 23andMe, the map argument is interpreted as a version number which we use to know the path
if args.is_23andMe:
	args.legend_fp = ("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/supplementary/23andMe_SNPs/v" + 
	                  str(args.legend_fp) + "/all_snp_info-" + 
	                  str(args.legend_fp) + ".txt")
	map_coord_name = "all.data.id"
	map_rs_name = "assay.name"
else:
	map_coord_name = "coord"
	map_rs_name = "rsid"

## READ DATA
# sumstat
sumstat = pd.read_csv(args.sumstat_fp, delim_whitespace=True)
	# if 23andMe is passed, drop all the missing rows so we only have the SNPs with effects
if args.is_23andMe:
	sumstat = sumstat.dropna(subset=['pvalue','effect'])
	sumstat = sumstat[sumstat.pass=="Y"]
# map
legend = pd.read_csv(args.legend_fp, delim_whitespace=True, 
                     usecols=[map_coord_name, map_rs_name])


## CLEAN AND MERGE DATA
# clean summary statistic coordinate column to ensure we get as many matches as possible
if not args.is_23andMe:
	# needed for parental longevity (deelen)
	sumstat[args.sumstat_coord_name] = sumstat[args.sumstat_coord_name].str.replace(":ID", "")
	# needed for osteoarthritis
	sumstat[args.sumstat_coord_name] = sumstat[args.sumstat_coord_name].str.replace("_.*$", "")

# merge summary statistics with legend
joined = pd.merge(sumstat, legend, 
                  left_on=args.sumstat_coord_name, right_on=map_coord_name, how="left")

sumstat_cols = sumstat.columns
del sumstat
del legend

joined = joined[sumstat_cols.append(pd.Index([map_rs_name]))]

joined = joined.rename({args.sumstat_coord_name: "unique_id", map_rs_name : "SNP"}, axis="columns")


## WRITE DATA
print("Warning: ", 
      sum(pd.isna(joined.SNP)), " of ", len(joined), 
      " SNP coordinates were not matched with RS numbers.")
joined.to_csv(args.output_fp, index=False, sep='\t')
