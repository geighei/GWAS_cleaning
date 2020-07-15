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
	legend_cols = [map_coord_name, map_rs_name, "alleles"]
else:
	map_coord_name = "coord"
	map_rs_name = "rsid"
	legend_cols = [map_coord_name, map_rs_name]

## READ DATA
# logging
print("Mapping RS numbers to summary statistic at file: " + args.sumstat_fp)
# sumstat
sumstat = pd.read_csv(args.sumstat_fp, delim_whitespace=True)
if args.is_23andMe:
	# if 23andMe is passed, drop all the missing rows so we only have the SNPs with effects that passed QC
	sumstat = sumstat.dropna(subset=['pvalue','effect'])
	sumstat = sumstat[sumstat["pass"]=="Y"]
# map
legend = pd.read_csv(args.legend_fp, delim_whitespace=True, usecols=legend_cols)

## CLEAN AND MERGE DATA
# clean summary statistic coordinate column to ensure we get as many matches as possible
if not args.is_23andMe:
	# needed for parental longevity (deelen)
	sumstat[args.sumstat_coord_name] = sumstat[args.sumstat_coord_name].str.replace(":ID", "")
	# needed for osteoarthritis
	sumstat[args.sumstat_coord_name] = sumstat[args.sumstat_coord_name].str.replace("_.*$", "")


# merge summary statistics with legend
print("Merging summary statistic with RS-number legend using columns " +
	args.sumstat_coord_name + " and " + map_coord_name + ".")
joined = pd.merge(sumstat, legend, 
                  left_on=args.sumstat_coord_name, right_on=map_coord_name, how="left")

sumstat_cols = sumstat.columns
del sumstat
del legend

# allele and snp cleaning for 23 and Me since the stats also don't have this info
if args.is_23andMe:
	# read in list of 1kG flipped SNPs that must be removed
	problematic_coords = pd.read_csv("/home/ubuntu/biroli/geighei/data/GWAS_sumstats/supplementary/23andMe_SNPs/1000G_flipped_all_SNPID.txt",
	                               header=None, delim_whitespace=True)
	# re-format chr and position columns so we have chr:pos column
	joined["chr_pos"] = joined["scaffold"].str.strip("chr") + ":" + legend["position"].astype(str))
	# identify set of coords in df that aren't problematic and filter data frame on them
	good_coords = set(joined.chr_pos).difference(problematic_coords[0])
	joined = joined[joined.chr_pos.isin(good_coords)]

	# "alleles" column is of form "A/C"; effect corresponds to reference allele which is "alphabetically greater" 
	# split and change legend_cols so we keep new allele cols
	joined[["A2","A1"]] = joined["alleles"].str.split('/',expand=True)
	legend_cols.remove("alleles")
	legend_cols.extend(["A1", "A2"])

# select only the columns we need and rename to standardize
legend_cols.remove(map_coord_name)
joined = joined[sumstat_cols.append(pd.Index(legend_cols))]
joined = joined.rename({args.sumstat_coord_name: "unique_id", map_rs_name : "SNP"}, axis="columns")


## WRITE DATA
print("Mapping of file completed. ", 
      sum(pd.isna(joined["SNP"]) | joined["SNP"].isnull()), 
      " of ", len(joined), " SNP coordinates were not matched with RS numbers.")
print("Writing output to " + args.output_fp + "." + "\n")
joined.to_csv(args.output_fp, index=False, sep='\t', na_rep="NA")
