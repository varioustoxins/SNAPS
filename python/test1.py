#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 10:45:35 2018

@author: aph516
"""

import numpy as np
import pandas as pd
from plotnine import *
from scipy.stats import norm
from scipy.optimize import linear_sum_assignment
from math import isnan, log10

#"~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt"

def import_obs_shifts(filename):
    #### Import the observed chemical shifts
    obs_long = pd.read_table(filename)
    obs_long = obs_long[["Residue_PDB_seq_code","Residue_label","Atom_name","Chem_shift_value"]]
    obs_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
    obs_long["SS_name"] = obs_long["Res_N"].astype(str) + obs_long["Res_type"]  # This needs to convert Res_type to single letter first
    obs_long = obs_long.reindex(columns=["Res_N","Res_type","SS_name","Atom_type","Shift"])
    
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
    tmp = obs_long[["Res_N","Res_type","SS_name"]]
    tmp = tmp.drop_duplicates(subset="SS_name")
    tmp.index = tmp["Res_N"]
    obs = pd.concat([tmp, obs], axis=1)
    
    obs.index = obs["SS_name"]
    
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
    
log_prob_matrix = calc_log_prob_matrix(obs, preds, sf=2)

def find_best_assignment(obs, preds, log_prob_matrix):
    # Use the Hungarian algorithm to find the highest probability matching 
    # (ie. the one with the lowest log probability sum)
    # Return a data frame with matched observed and predicted shifts, and the raw matching
    valid_atoms = ["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]
    
    row_ind, col_ind = linear_sum_assignment(log_prob_matrix)
    
    obs_names = [log_prob_matrix.index[r] for r in row_ind]
    pred_names = [log_prob_matrix.columns[c] for c in col_ind]
    #return(list(zip(obs_names, pred_names)))
    
    assign_df = pd.DataFrame({
            "Res_name":pred_names,
            "SS_name":obs_names
            })
    #assign_df.index = assign_df["Res_name"]
    
    assign_df = pd.merge(assign_df, preds[["Res_N","Res_type"]], on="Res_name")
    assign_df = assign_df[["Res_name","Res_N","Res_type","SS_name"]]
    assign_df = pd.merge(assign_df, obs.loc[:, obs.columns.isin(valid_atoms+["SS_name"])], on="SS_name")
    # Above line raises an error about index/column confusion, which needs fixing.
    assign_df = pd.merge(assign_df, preds.loc[:, preds.columns.isin(valid_atoms)], on="Res_name", suffixes=("","_pred"))
    
    assign_df = assign_df.sort_values(by="Res_N")
    
    return(assign_df, [row_ind, col_ind])

assign_df, matching = find_best_assignment(obs, preds, log_prob_matrix)

def plot_strips(assign_df, atom_list=["C","Cm1","CA","CAm1","CB","CBm1"]):
    # Make a strip plot of the assignment, using only the atoms in atom_list
    
    # First, convert assign_df from wide to long
    plot_df = assign_df.loc[:,["Res_N", "Res_type", "Res_name", "SS_name"]+atom_list]
    plot_df = plot_df.melt(id_vars=["Res_N", "Res_type", "Res_name", "SS_name"],
                               value_vars=atom_list, var_name="Atom_type", value_name="Shift")
    
    # Add columns with information to be plotted
    plot_df["i"] = "0"     # This column determines if shift is from the i or i-1 residue
    plot_df.loc[plot_df["Atom_type"].isin(["Cm1","CAm1","CBm1"]),"i"] = "-1"
    plot_df["Atom_type"] = plot_df["Atom_type"].replace({"Cm1":"C", "CAm1":"CA", "CBm1":"CB"}) # Simplify atom type
    
    plot_df["seq_group"] = plot_df["Res_N"] + plot_df["i"].astype("int")
    
    # Pad Res_name column with spaces so that sorting works correctly
    plot_df["Res_name"] = plot_df["Res_name"].str.pad(6)
    plot_df["x_name"] = plot_df["Res_name"] + "_(" + plot_df["SS_name"] + ")"
    
    # Make the plot
    plt = ggplot(aes(x="x_name"), data=plot_df) + geom_point(aes(y="Shift", colour="i"))
    plt = plt + geom_line(aes(y="Shift", group="seq_group"))        # Add lines connecting i to i-1 points
    plt = plt + geom_line(aes(y="Shift", group="Res_N"), linetype="dashed")
    plt = plt + facet_grid("Atom_type~.", scales="free") + scale_colour_brewer(type="Qualitative", palette="Set1") 
    plt = plt + xlab("Residue name") + ylab("Chemical shift (ppm)")
    plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90)) + scale_y_reverse()
    
    return(plt)
    
plot_strips(assign_df.iloc[0:10, :])
