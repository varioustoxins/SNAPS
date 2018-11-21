#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:30:40 2018

@author: aph516
"""

import pandas as pd
from subprocess import run

path = "/Users/aph516/GitHub/NAPS/"


#%%
#### Test all proteins in the testset

testset_df = pd.read_table(path+"data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = path+"data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"
testset_df["preds_file"] = path+"data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"
testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
testset_df.index = testset_df["ID"]

for i in testset_df.index:
    print(testset_df.loc[i, "obs_file"])
    args = [path+"python/NAPS.py", testset_df.loc[i, "obs_file"], testset_df.loc[i, "preds_file"], 
            path+"output/testset/"+testset_df.loc[i, "out_name"]+".txt"]
    run(args)



#%%
#### Test alt assignments