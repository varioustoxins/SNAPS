#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 10:45:35 2018

@author: aph516
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from math import isnan, log10

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
    obs = pd.merge(obs, obs_m1, how="left", left_index=True, right_index=True)
    
    # Restrict to specific atom types
    atom_list = ["H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"]
    obs = obs[atom_list]
    
    # Add the other data back in
    tmp = obs_long[["Res_N","Res_type","Res_name"]]
    tmp = tmp.drop_duplicates(subset="Res_name")
    tmp.index = tmp["Res_N"]
    obs = pd.concat([tmp, obs], axis=1)
    
    obs.index = obs["Res_name"]
    
    return(obs)

obs = import_obs_shifts("~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt")

def read_shiftx2(input_file, offset=0):
    #### Import the predicted chemical shifts
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
    
    preds.index = preds["Res_name"]
    
    return(preds)

preds = read_shiftx2("~/GitHub/NAPS/data/testset/shiftx2_results/A001_1KF3A.cs")

#### Create a probability matrix
obs1=obs.iloc[0]
pred1=preds.iloc[0]

def calc_match_probability(obs1, pred1,
                           atom_set=set(["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]), 
                           atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231, 
                                    'C':0.5330, 'CA':0.4412, 'CB':0.5163, 
                                    'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}, sf=1, default_prob=0.01):
    # Calculate a match score between an observed spin system and a predicted residue, assuming a Gaussian probability distribution.
    # obs1 and obs2 are single rows from the obs and preds data frames. Note that they become type Series
    # default_prob is the probability assigned when an observation or prediction is missing
    # atom_set is a set used to restrict to only certain measurements
    # atom_sd is the expected standard deviation for each atom type
    # sf is a scaling factor for the entire atom_sd dictionary
    
    # At the moment, this doesn't deal with case where both obs and pres are missing an atom.
    
     # If observed residue has an amide proton, it can't be proline.
    if pred1["Res_type"]=="P" and not isnan(obs1["H"]):      # I don't think this works properly
        return(0)
    else:
        # Throw away any non-atom columns
        obs1 = obs1.loc[atom_set.intersection(obs1.index)]
        pred1 = pred1.loc[atom_set.intersection(pred1.index)]
        df = pd.DataFrame({'obs':obs1, 'pred':pred1})   # Make obs1 and pred1 into columns of a data frame
        df["Delta"] = df['obs'] - df['pred']            # Calculate difference between observation and prediction  
        df["Prob"] = 1
        for i in df.index:          # For each atom type present, calculate the probability that Delta would be this size or larger
            if isnan(df.loc[i, "Delta"]):   # Use default probability if data is missing.
                df.loc[i, "Prob"] = default_prob
            else:       # Otherwise, probability taken from cumulative distribution function.
                df.loc[i, "Prob"] = 2*norm.cdf(-abs(df.loc[i, "Delta"]), scale=atom_sd[i]*sf)
        
        overall_prob = df["Prob"].prod(skipna=False)
        return(overall_prob)
    
calc_match_probability(obs.loc["1LYS"], preds.loc["1K"], sf=2)

def calc_probability_matrix(obs, pred,
                            atom_set=set(["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]), 
                            atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231, 
                                    'C':0.5330, 'CA':0.4412, 'CB':0.5163, 
                                    'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}, sf=1, default_prob=0.01):
    # Calculate a matrix of match probabilities
    prob_matrix = pd.DataFrame(np.NaN, index=obs.index, columns=preds.index)    # Initialise matrix as NaN
    
    for i in obs.index:
        print(i)
        for j in preds.index:
            prob_matrix.loc[i, j] = calc_match_probability(obs.loc[i,:], preds.loc[j,:], atom_set, atom_sd, sf, default_prob)
            
    return(prob_matrix)
    
#prob_matrix = calc_probability_matrix(obs, preds)

def calc_log_prob_matrix(obs, pred,
                            atom_set=set(["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]), 
                            atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231, 
                                    'C':0.5330, 'CA':0.4412, 'CB':0.5163, 
                                    'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}, sf=1, default_prob=0.01):
    # Calculate a matrix of -log10(match probabilities)
    prob_matrix = pd.DataFrame(np.NaN, index=obs.index, columns=preds.index)    # Initialise matrix as NaN
    
    for i in obs.index:
        print(i)
        for j in preds.index:
            prob_matrix.loc[i, j] = calc_match_probability(obs.loc[i,:], preds.loc[j,:], atom_set, atom_sd, sf, default_prob)
    
    # Calculate log of matrix
    log_prob_matrix = -prob_matrix[prob_matrix>0].applymap(log10)
    log_prob_matrix[log_prob_matrix.isna()] = 2*np.nanmax(log_prob_matrix.values)
        
    return(log_prob_matrix)
    
log_prob_matrix = calc_log_prob_matrix(obs, preds)
