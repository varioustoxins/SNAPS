#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 14:30:29 2019

@author: aph516
"""


import numpy as np
import pandas as pd
from SNAPS_importer import SNAPS_importer
from SNAPS_assigner import SNAPS_assigner
from pathlib import Path
from scipy.stats import norm
from copy import deepcopy
from math import isnan, log10

path = Path("..")
#path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/Alex/GitHub/NAPS")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")

testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/"data/testset/simplified_BMRB"/file for file in testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/"data/testset/shiftx2_results"/file for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df["noshifty_file"] = [path/"data/testset/noshifty_results"/file for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df.index = testset_df["ID"]

#%%

# Import observed and predicted shifts
id = "A003"

importer = SNAPS_importer()
importer.import_testset_shifts(testset_df.loc[id, "obs_file"])


#%%
a = SNAPS_assigner()

a.obs = importer.obs

# Import config file
a.read_config_file(path/"config/config_yaml.txt")

tmp = a.import_pred_shifts(testset_df.loc[id, "preds_file"], "shiftx2")

# Do the analysis
tmp = a.prepare_obs_preds()
tmp = a.calc_log_prob_matrix()
tmp = a.calc_mismatch_matrix()
tmp = a.assign_from_preds()
tmp = a.add_consistency_info(threshold=a.pars["seq_link_threshold"])

tmp = a.find_consistent_assignments()

obs = a.obs
preds = a.preds
log_prob_matrix = a.log_prob_matrix
mismatch_matrix = a.mismatch_matrix
assign_df = a.assign_df
