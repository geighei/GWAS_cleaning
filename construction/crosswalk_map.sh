#!/usr/bin/env bash
shopt -s extglob
##
## Usage: 
# ./crosswalk_map.sh --fp ~/biroli/geighei/data/GWAS_sumstats/construction/height/height_pheno \
# 	--cross ~/biroli/geighei/data/GWAS_sumstats/construction/test.crosswalk


while [[ "$#" -gt 0 ]]; do
    case $1 in
    	--fp) fp="$2"; shift ;;
        --cross) crosswalk="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

awk 'BEGIN {print "FID", "IID", "resid"}
NR==FNR {if($2 != "NA") {a[$1] = $2 + 1; next}}
NR!=FNR {if(a[$2] != 0) {print $1, $1, a[$2] - 1}
else if(a[$2] == 0) {print $1, $1, "NA"}}' \
	${fp}.txt ${crosswalk} > ${fp}.PREPARED.pheno