#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 10:45:35 2018

@author: aph516
"""

import numpy as np
import pandas as pd

#"~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt"

def import_obs_shifts(filename):
    #### Import the observed chemical shifts
    obs_long = pd.read_table(filename)
    obs_long = obs_long[["Residue_PDB_seq_code","Residue_label","Atom_name","Chem_shift_value"]]
    obs_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
    obs_long["Res_name"] = obs_long["Res_N"].astype(str) + obs_long["Res_type"]  # This needs to convert Res_type to single letter first
    obs_long = obs_long.reindex(columns=["Res_N","Res_type","Res_name","Atom_type","Shift"])
    
    # Convert from long to wide
    obs = obs_long.pivot(index="Res_N", columns="Atom_type", values="Shift")
    
    # Make columns for the i-1 observed shifts of C, CA and CB
    obs_m1 = obs[["C","CA","CB"]]
    obs_m1.index = obs_m1.index+1
    obs_m1.columns = ["Cm1","CAm1","CBm1"]
    obs = pd.concat([obs, obs_m1], axis=1)
    
    # Restrict to specific atom types
    atom_list = ["H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"]
    obs = obs[atom_list]
    
    return(obs)

obs2 = import_obs_shifts("~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt")

#### Import the predicted chemical shifts
preds = pd.read_csv("~/GitHub/NAPS/data/testset/shiftx2_results/A001_1KF3A.cs")
preds["Res_name"] = preds["NUM"].astype(str)+preds["RES"]
if any(preds.columns == "CHAIN"):   preds = preds.drop("CHAIN", axis=1)     # Assuming that there's only one CHAIN in the predictions...
preds = preds.reindex(columns=["NUM","RES","Res_name","ATOMNAME","SHIFT"])  
preds.columns = ["Res_N","Res_type","Res_name","Atom_type","Shift"]

# Convert from wide to long format
preds_wide = preds.pivot(index="Res_N", columns="Atom_type", values="Shift")

# Make columns for the i-1 predicted shifts of C, CA and CB
preds_wide_m1 = preds_wide[["C","CA","CB"]].copy()
preds_wide_m1.index = preds_wide_m1.index+1
#preds_wide_m1 = preds_wide_m1[["C","CA","CB"]]
preds_wide_m1.columns = ["Cm1","CAm1","CBm1"]
preds_wide = pd.merge(preds_wide, preds_wide_m1, how="left", left_index=True, right_index=True)

# Restrict to only certain atom types
atom_list = ["H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"]
preds_wide = preds_wide[atom_list]

# Add the other data back in
tmp = preds[["Res_N","Res_type","Res_name"]]
tmp = tmp.drop_duplicates(subset="Res_name")
tmp.index = tmp["Res_N"]
preds_wide = pd.concat([tmp, preds_wide], axis=1)

def read_shiftx2(input_file, offset=0):
    preds_long = pd.read_csv(input_file)
    preds_long["NUM"] = preds_long["NUM"] + offset    # Apply any offset to residue numbering
    preds_long["Res_name"] = preds_long["NUM"].astype(str)+preds_long["RES"]
    if any(preds_long.columns == "CHAIN"):   preds_long = preds_long.drop("CHAIN", axis=1)     # Assuming that there's only one CHAIN in the predictions...
    preds_long = preds_long.reindex(columns=["NUM","RES","Res_name","ATOMNAME","SHIFT"])  
    preds_long.columns = ["Res_N","Res_type","Res_name","Atom_type","Shift"]
    
    # Convert from wide to long format
    preds = preds_long.pivot(index="Res_N", columns="Atom_type", values="Shift")
    
    # Make columns for the i-1 predicted shifts of C, CA and CB
    preds_m1 = preds[["C","CA","CB"]].copy()
    preds_m1.index = preds_m1.index+1
    #preds_m1 = preds_m1[["C","CA","CB"]]
    preds_m1.columns = ["Cm1","CAm1","CBm1"]
    preds = pd.merge(preds, preds_m1, how="left", left_index=True, right_index=True)
    
    # Restrict to only certain atom types
    atom_list = ["H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"]
    preds = preds[atom_list]
    
    # Add the other data back in
    tmp = preds_long[["Res_N","Res_type","Res_name"]]
    tmp = tmp.drop_duplicates(subset="Res_name")
    tmp.index = tmp["Res_N"]
    preds = pd.concat([tmp, preds], axis=1)
    
    return(preds)

preds = read_shiftx2("~/GitHub/NAPS/data/testset/shiftx2_results/A001_1KF3A.cs")
