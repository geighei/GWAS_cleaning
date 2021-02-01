import pandas as pd
import numpy as np

map_1kg = pd.read_csv("1kGenomes/legends/1kGenomes_allCHR.legend", delim_whitespace=True)
map_hrc = pd.read_csv("HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map", delim_whitespace=True)
map_hrc = map_hrc.rename({"chr":"CHR", "pos":"POS"}, axis="columns")
map_hrc["coord"] = map_hrc.CHR.map(str) + ":" + map_hrc.POS.map(str)

rsid_map = pd.concat([map_1kg, map_hrc], ignore_index=True, sort=False)
rsid_map.drop_duplicates(subset=['rsid'], inplace=True)
rsid_map.to_csv("rsid_map.txt", index=False, sep='\t')
