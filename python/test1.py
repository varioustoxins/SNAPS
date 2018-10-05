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
from Bio.SeqUtils import seq1

#"~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt"

def import_obs_shifts(filename, remove_Pro=True):
    #### Import the observed chemical shifts
    obs_long = pd.read_table(filename)
    obs_long = obs_long[["Residue_PDB_seq_code","Residue_label","Atom_name","Chem_shift_value"]]
    obs_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
    obs_long["Res_type"] = obs_long["Res_type"].apply(seq1)
    obs_long["SS_name"] = obs_long["Res_N"].astype(str) + obs_long["Res_type"]  # This needs to convert Res_type to single letter first
    obs_long = obs_long.reindex(columns=["Res_N","Res_type","SS_name","Atom_type","Shift"])
    
    # Convert from long to wide
    obs = obs_long.pivot(index="Res_N", columns="Atom_type", values="Shift")
    
    # Make columns for the i-1 observed shifts of C, CA and CB
    obs_m1 = obs[list({"C","CA","CB"}.intersection(obs.columns))]
    obs_m1.index = obs_m1.index+1
    obs_m1.columns = obs_m1.columns + "m1"
    obs = pd.merge(obs, obs_m1, how="left", left_index=True, right_index=True)
    
    # Restrict to specific atom types
    atom_set = {"H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"}
    obs = obs[list(atom_set.intersection(obs.columns))]
    
    # Add the other data back in
    tmp = obs_long[["Res_N","Res_type","SS_name"]]
    tmp = tmp.drop_duplicates(subset="SS_name")
    tmp.index = tmp["Res_N"]
    obs = pd.concat([tmp, obs], axis=1)
    
    obs.index = obs["SS_name"]
    
    if remove_Pro:  obs = obs.drop(obs.index[obs["Res_type"]=="PRO"]) # Remove prolines, as they wouldn't be observed in a real spectrum
    
    return(obs)


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
    preds_m1 = preds[list({"C","CA","CB"}.intersection(preds.columns))].copy()
    preds_m1.index = preds_m1.index+1
    preds_m1.columns = preds_m1.columns + "m1"
    preds = pd.merge(preds, preds_m1, how="left", left_index=True, right_index=True)
    # TODO: also do this for Res_type
    
    # Restrict to only certain atom types
    atom_set = {"H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"}
    preds = preds[list(atom_set.intersection(preds.columns))]
    
    # Add the other data back in
    tmp = preds_long[["Res_N","Res_type","Res_name"]]
    tmp = tmp.drop_duplicates(subset="Res_name")
    tmp.index = tmp["Res_N"]
    preds = pd.concat([tmp, preds], axis=1)
    
    preds.index = preds["Res_name"]
    
    return(preds)


def add_dummy_rows(obs, preds):
    # Add dummy rows to obs and preds to bring them to the same length
    # Also throw away any atom types that aren't present in both obs and preds
    
    # Delete any prolines in preds
    preds = preds.drop(preds.index[preds["Res_type"]=="P"])
    
    # Restrict atom types
    atom_set = {"H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"}
    obs_metadata = list(set(obs.columns).difference(atom_set))
    preds_metadata = list(set(preds.columns).difference(atom_set))
    shared_atoms = list(atom_set.intersection(obs.columns).intersection(preds.columns))
    obs = obs[obs_metadata+shared_atoms]
    preds = preds[preds_metadata+shared_atoms]
    
    # Create columns to keep track of dummy status
    preds["Dummy_res"] = False
    obs["Dummy_SS"] = False
    
    N = len(obs.index)
    M = len(preds.index)
    
    if N>M:     # If there are more spin systems than predictions
        dummies = pd.DataFrame(np.NaN, index=["dummy_res_"+str(i) for i in 1+np.arange(N-M)], columns = preds.columns)
        dummies["Res_name"] = dummies.index
        dummies["Dummy_res"] = True
        preds = preds.append(dummies)        
    elif M>N:
        dummies = pd.DataFrame(np.NaN, index=["dummy_SS_"+str(i) for i in 1+np.arange(M-N)], columns = obs.columns)
        dummies["SS_name"] = dummies.index
        dummies["Dummy_SS"] = True
        obs = obs.append(dummies)
        #obs.loc[["dummy_"+str(i) for i in 1+np.arange(M-N)]] = np.NaN
        #obs.loc[obs.index[N:M], "SS_name"] = ["dummy_"+str(i) for i in 1+np.arange(M-N)]
    
    return(obs, preds)


def calc_match_probability(obs, pred1,
                           atom_set=set(["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]), 
                           atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231, 
                                    'C':0.5330, 'CA':0.4412, 'CB':0.5163, 
                                    'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}, sf=1, default_prob=0.01):
    # Calculate match scores between all observed spin systems and a single predicted residue
    # default_prob is the probability assigned when an observation or prediction is missing
    # atom_set is a set used to restrict to only certain measurements
    # atom_sd is the expected standard deviation for each atom type
    # sf is a scaling factor for the entire atom_sd dictionary
    
    # Doesn't currently deal specially with prolines 
    
    # Throw away any non-atom columns
    obs = obs.loc[:, atom_set.intersection(obs.columns)]
    pred1 = pred1.loc[atom_set.intersection(pred1.index)]
    
    # Calculate shift differences and probabilities for each observed spin system
    delta = obs - pred1
    prob = delta.copy()
    prob.iloc[:,:] = 1
    
    # Make a note of NA positions in delta, and set them to zero (this avoids warnings when using norm.cdf later)
    na_mask = delta.isna()
    delta[na_mask] = 0
    for c in delta.columns:
        # Use the cdf to calculate the probability of a delta *at least* as great as the actual one
        prob[c] = 2*norm.cdf(-1*abs(pd.to_numeric(delta[c])), scale=atom_sd[c]*sf)
    
    # In positions where data was missing, use a default probability
    prob[na_mask] = default_prob
    
    # Calculate overall probability of each row
    overall_prob = prob.prod(skipna=False, axis=1)
    return(overall_prob)

#calc_match_probability(obs, pred1)

def calc_log_prob_matrix(obs, preds,
                            atom_set=set(["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]), 
                            atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231, 
                                    'C':0.5330, 'CA':0.4412, 'CB':0.5163, 
                                    'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}, sf=1, default_prob=0.01, verbose=False):
    # Calculate a matrix of -log10(match probabilities)
    prob_matrix = pd.DataFrame(np.NaN, index=obs.index, columns=preds.index)    # Initialise matrix as NaN
    
    for i in preds.index:
        if verbose: print(i)
        prob_matrix.loc[:, i] = calc_match_probability(obs, preds.loc[i,:], atom_set, atom_sd, sf, default_prob)
    
    # Calculate log of matrix
    log_prob_matrix = -prob_matrix[prob_matrix>0].applymap(log10)
    log_prob_matrix[log_prob_matrix.isna()] = 2*np.nanmax(log_prob_matrix.values)
    log_prob_matrix.loc[obs["Dummy_SS"], :] = 0
    log_prob_matrix.loc[:, preds["Dummy_res"]] = 0
        
    return(log_prob_matrix)


def find_best_assignment(obs, preds, log_prob_matrix):
    # Use the Hungarian algorithm to find the highest probability matching 
    # (ie. the one with the lowest log probability sum)
    # Return a data frame with matched observed and predicted shifts, and the raw matching
    valid_atoms = ["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"]
    
    row_ind, col_ind = linear_sum_assignment(log_prob_matrix)
    
    obs.index.name = "SS_index"
    preds.index.name = "Res_index"
    
    #Create assignment dataframe
    obs_names = [log_prob_matrix.index[r] for r in row_ind]
    pred_names = [log_prob_matrix.columns[c] for c in col_ind]
    
    assign_df = pd.DataFrame({
            "Res_name":pred_names,
            "SS_name":obs_names
            })
    
    # Merge residue information, shifts and predicted shifts into assignment dataframe
    assign_df = pd.merge(assign_df, preds[["Res_N","Res_type", "Res_name", "Dummy_res"]], on="Res_name")
    assign_df = assign_df[["Res_name","Res_N","Res_type","SS_name", "Dummy_res"]]
    assign_df = pd.merge(assign_df, obs.loc[:, obs.columns.isin(valid_atoms+["SS_name","Dummy_SS"])], on="SS_name")
        # Above line raises an error about index/column confusion, which needs fixing.
    assign_df = pd.merge(assign_df, preds.loc[:, preds.columns.isin(valid_atoms+["Res_name"])], on="Res_name", suffixes=("","_pred"))
    
    assign_df = assign_df.sort_values(by="Res_N")
    assign_df["Log_prob"] = log_prob_matrix.lookup(log_prob_matrix.index[row_ind], log_prob_matrix.columns[col_ind])
    
    return(assign_df, [row_ind, col_ind])


def find_alt_assignments(log_prob_matrix, best_match_indexes, by_res=True,  verbose=False):
    # Function to find the second best assignment for each residue or spin system
    # Does this by setting the log probability to a very high value for each residue in turn, and rerunning the assignment
    # best_match_indexes is the [row_ind, col_ind] output from find_best_assignment()
    # If by_res is true, calculate next best assignment for each residue. Otherwise, calculate it for each spin system.
    # Outputs a dataframe with the next best matching for each residue/spin system, and a list of their relative probabilities
    
    # Calculate sum probability for the best matching
    best_sum_prob = sum(log_prob_matrix.lookup(log_prob_matrix.index[best_match_indexes[0]], 
                               log_prob_matrix.columns[best_match_indexes[1]]))
    
    # Convert best_match_indexes into a list of (obs.index, preds.index) tuples
    best_matching = list(zip(obs.index[best_match_indexes[0]], preds.index[best_match_indexes[1]]))
    
    if by_res:
        # Create a dataframe to store the results
        matching_df = pd.DataFrame(index=obs.index[best_match_indexes[0]], columns=preds.index)
        rel_log_prob = pd.Series(index = preds.index) 
        penalty = 2*log_prob_matrix.max().max()     # This is the value used to penalise the best match for each residue
        
        for i in best_matching:
            if verbose: print(i)
            ss, res = i
            tmp = log_prob_matrix.copy()
            tmp.loc[ss, res] = penalty
            row_ind, col_ind = linear_sum_assignment(tmp)
            matching_df.loc[:,res] = preds.index[col_ind]
            rel_log_prob[res] = sum(tmp.lookup(tmp.index[row_ind], tmp.columns[col_ind])) - best_sum_prob
    
    return(matching_df, rel_log_prob)



    def check_assignment_consistency(assign_df, threshold=0.1):
        # Add columns to the assignment dataframe with the maximum mismatch to a sequential residue, and number of 'significant' mismatches
        # Any sequential residues with a difference greater than threshold count as mismatched
        
        # First check if there are any sequential atoms
        carbons = pd.Series(["C","CA","CB"])
        carbons_m1 = carbons + "m1"
        seq_atoms = carbons[carbons.isin(assign_df.columns) & carbons_m1.isin(assign_df.columns)]
        seq_atoms_m1 = seq_atoms+"m1"
        #seq_atoms = list(seq_atoms)
    
        if seq_atoms.size==0:
            # You can't do a comparison
            assign_df[["Max_mismatch_prev", "Max_mismatch_next", "Num_good_links_prev", "Num_good_links_next"]] = np.NaN
            return(assign_df)
        else:
            # First, get the i and i-1 shifts for the preceeding and succeeding residues
            tmp = assign_df.copy()
            tmp = tmp.loc[tmp["Dummy_res"]==False,]
            tmp.index = tmp["Res_N"]
            tmp = tmp[list(seq_atoms)+list(seq_atoms_m1)]
            tmp_next = tmp.copy()
            tmp_next.index -= 1
            tmp_prev = tmp.copy()
            tmp_prev.index += 1
            tmp = tmp.join(tmp_next, rsuffix="_next")
            tmp = tmp.join(tmp_prev, rsuffix="_prev")
            # Calculate mismatch for each atom type
            for atom in seq_atoms:
                tmp["d"+atom+"_prev"] = tmp[atom+"m1"] - tmp[atom+"_prev"]
                tmp["d"+atom+"_next"] = tmp[atom] - tmp[atom+"m1_next"]
            # Calculate maximum mismatch
            tmp["Max_mismatch_prev"] = tmp["d"+seq_atoms+"_prev"].max(axis=1, skipna=True)
            tmp["Max_mismatch_next"] = tmp["d"+seq_atoms+"_next"].max(axis=1, skipna=True)
            
            # Calculate number of consistent matches
            #
            tmp["Num_good_links_prev"] = (tmp["d"+seq_atoms+"_prev"]<threshold).sum(axis=1)
            tmp["Num_good_links_next"] = (tmp["d"+seq_atoms+"_next"]<threshold).sum(axis=1)
            #tmp["Num_good_links"] = tmp["Num_good_links_prev"] + tmp["Num_good_links_next"]
            
            # Join relevant columns back onto assign_df
            tmp["Res_N"] = tmp.index
            assign_df = assign_df.join(tmp.loc[:,["Max_mismatch_prev", "Max_mismatch_next", "Num_good_links_prev", "Num_good_links_next"]], on="Res_N")
           
            return(assign_df)
    

def plot_strips(assign_df, atom_list=["C","Cm1","CA","CAm1","CB","CBm1"]):
    # Make a strip plot of the assignment, using only the atoms in atom_list
    
    # Narrow down atom list to those actually present
    atom_list = list(set(atom_list).intersection(assign_df.columns))
    
    # First, convert assign_df from wide to long
    plot_df = assign_df.loc[:,["Res_N", "Res_type", "Res_name", "SS_name", "Dummy_res", "Dummy_SS"]+atom_list]
    plot_df = plot_df.melt(id_vars=["Res_N", "Res_type", "Res_name", "SS_name", "Dummy_res", "Dummy_SS"],
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
    plt = ggplot(aes(x="x_name"), data=plot_df) + geom_point(aes(y="Shift", colour="i", shape="Dummy_res"))
    plt = plt + scale_y_reverse() + scale_shape_manual(values=["o","x"])
    plt = plt + geom_line(aes(y="Shift", group="seq_group"), data=plot_df.loc[~plot_df["Dummy_res"],])        # Add lines connecting i to i-1 points
    plt = plt + geom_line(aes(y="Shift", group="x_name"), linetype="dashed")
    plt = plt + facet_grid("Atom_type~.", scales="free") + scale_colour_brewer(type="Qualitative", palette="Set1") 
    plt = plt + xlab("Residue name") + ylab("Chemical shift (ppm)")
    plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90))
    
    return(plt)

def plot_seq_mismatch(assign_df):
    # Make a plot of the maximum sequential mismatch between i-1,i and i+1 residues
    
    # Check that the assignment data frame has the right columns
    if not all(pd.Series(['Max_mismatch_prev', 'Max_mismatch_next']).isin(assign_df.columns)):
        return(None)
    else:
        # Pad Res_name column with spaces so that sorting works correctly
        assign_df["Res_name"] = assign_df["Res_name"].str.pad(6)
        assign_df["x_name"] = assign_df["Res_name"] + "_(" + assign_df["SS_name"] + ")"
        
        # Make the plot
        plt = ggplot(aes(x="x_name"), data=assign_df) 
        plt = plt + geom_col(aes(y="abs(Max_mismatch_prev)"))
        plt = plt + xlab("Residue name") + ylab("Mismatch to previous residue (ppm)")
        plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90))
               
        return(plt)

def NAPS_single(obs_file, preds_file, out_name, 
                out_path="../output", plot_path="../plots", 
                use_atoms=["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"], make_plots=True):
    # Function to assign a single protein from files of observed and predicted shifts
    
    # Import the observed and predicted shifts
    obs = import_obs_shifts(obs_file)
    preds = read_shiftx2(preds_file)
    
    # Add dummy rows so that obs and preds are the same length
    obs, preds = add_dummy_rows(obs, preds)
    
    # Calculate the log probability for each observation/prediction pair
    log_prob_matrix = calc_log_prob_matrix(obs, preds, sf=2)
    
    # Find the assignment with the highest overall probability
    assign_df, matching = find_best_assignment(obs, preds, log_prob_matrix)
    assign_df = check_assignment_consistency(assign_df, threshold=0.05)
    assign_df.to_csv(out_path+"/"+out_name+".txt", sep="\t", na_rep="NA", float_format="%.3f")
    
    # Produce a strip plot
    if make_plots:
        plt = plot_strips(assign_df.iloc[:, :])
        plt.save(out_name+".pdf", path=plot_path, height=210, width=max(297, preds.index.size*3), units="mm", limitsize=False)
    return(1)
    
def NAPS_batch(file_df, out_path="../output", plot_path="../plots",
               use_atoms=["H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"], make_plots=True):
    # Assign multiple proteins.
    # file_df is a pandas DataFrame containing obs_file, preds_file and out_name columns
    for i in file_df.index:
        print(file_df.loc[i, "ID"])
        NAPS_single(file_df.loc[i, "obs_file"], file_df.loc[i, "preds_file"], 
                    file_df.loc[i, "out_name"], out_path, plot_path, use_atoms, make_plots)
    return(1)

#### Main
    
#NAPS_single("~/GitHub/NAPS/data/testset/simplified_BMRB/6338.txt", "~/GitHub/NAPS/data/testset/shiftx2_results/A002_1XMTA.cs", "A002_6338")    

# Assign a batch of proteins
testset_df = pd.read_table("../data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = "../data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"
testset_df["preds_file"] = "../data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"
testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
testset_df.index = testset_df["ID"]

#NAPS_batch(testset_df.iloc[:,:], "../output/testset", "../plots/testset", make_plots=True)

i = "A002"
obs = import_obs_shifts("~/GitHub/NAPS/data/testset/simplified_BMRB/"+testset_df.BMRB[i].astype(str)+".txt")
preds = read_shiftx2("~/GitHub/NAPS/data/testset/shiftx2_results/"+i+"_"+testset_df.PDB[i]+".cs")
obs, preds = add_dummy_rows(obs, preds)

#### Create a probability matrix
#obs1=obs.iloc[0]
#pred1=preds.iloc[0]
log_prob_matrix = calc_log_prob_matrix(obs, preds, sf=2, verbose=False)
assign_df, best_match_indexes = find_best_assignment(obs, preds, log_prob_matrix)
assign_df = check_assignment_consistency(assign_df)
alt_assignments, alt_assignment_scores = find_alt_assignments(log_prob_matrix, best_match_indexes)
#assign_df.to_csv("../output/A001_4032.txt", sep="\t", float_format="%.3f")

plot_strips(assign_df.iloc[:, :])
#plt = plot_seq_mismatch(assign_df)
#plt.save("A001_4032_mismatch.pdf", path="../plots", height=210, width=297, units="mm")
