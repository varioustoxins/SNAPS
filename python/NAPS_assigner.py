#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:36:07 2018

@author: aph516
"""

import numpy as np
import pandas as pd
from plotnine import *
from scipy.stats import norm, multivariate_normal
from scipy.optimize import linear_sum_assignment
from math import isnan, log10, sqrt
from copy import deepcopy
from Bio.SeqUtils import seq1
import logging

class NAPS_assigner:
    # Attributes
    obs = None
    preds = None
    log_prob_matrix = None
    assign_df = None
    alt_assign_df = None
    best_match_indexes = None
    pars = {"pred_offset": 0,
            "atom_set": {"H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"},
            "atom_sd": {'H':0.1711, 'N':1.1169, 'HA':0.1231,
                        'C':0.5330, 'CA':0.4412, 'CB':0.5163,
                        'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}}
    
    # Functions
    def import_pred_shifts(self, input_file, filetype, offset=None):
        """ Import predicted chemical shifts from a ShiftX2 results file.
        
        filetype: either "shiftx2" or "sparta+"
        offset: an optional integer to add to the ShiftX2 residue number.
        """
        
        # If no offset value is defined, use the default one
        if offset==None:
            offset = self.pars["pred_offset"]
        
        if filetype == "shiftx2":
            preds_long = pd.read_csv(input_file)
            if any(preds_long.columns == "CHAIN"):
                if len(preds_long["CHAIN"].unique())>1:
                    print("Chain identifier dropped - if multiple chains are "+
                          "present in the predictions, they will be merged.")
                preds_long = preds_long.drop("CHAIN", axis=1)     
            preds_long = preds_long.reindex(columns=["NUM","RES","ATOMNAME",
                                                     "SHIFT"])  
            preds_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
        elif filetype == "sparta+":
            # Work out where the column names and data are
            with open(input_file, 'r') as f:
                for num, line in enumerate(f, 1):
                    if line.find("VARS")>-1:
                        colnames_line = num
                        colnames = line.split()[1:]
                        
            preds_long = pd.read_table(input_file, sep="\s+", names=colnames,
                                       skiprows=colnames_line+1)
            preds_long = preds_long[["RESID","RESNAME","ATOMNAME","SHIFT"]]
            preds_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
            
            # Sparta+ uses HN for backbone amide proton - convert to H
            preds_long.loc[preds_long["Atom_type"]=="HN", "Atom_type"] = "H"
        else:
            print("import_pred_shifts: invalid filetype '%s'." % (filetype))
            return(None)
        
        # Add sequence number offset and create residue names
        preds_long["Res_N"] = preds_long["Res_N"] + offset
        preds_long.insert(1, "Res_name", (preds_long["Res_N"].astype(str) + 
                  preds_long["Res_type"]))
        preds_long["Res_name"] = [s.rjust(5) for s in preds_long["Res_name"]]
            
        # Convert from wide to long format
        preds = preds_long.pivot(index="Res_N", columns="Atom_type", 
                                 values="Shift")
        
        # Add the other data back in
        tmp = preds_long[["Res_N","Res_type","Res_name"]]
        tmp = tmp.drop_duplicates(subset="Res_name")
        tmp.index = tmp["Res_N"]
        preds = pd.concat([tmp, preds], axis=1)
        
        # Make columns for the i-1 predicted shifts of C, CA and CB
        preds_m1 = preds[list({"C","CA","CB","Res_type"}.
                              intersection(preds.columns))].copy()
        preds_m1.index = preds_m1.index+1
        preds_m1.columns = preds_m1.columns + "m1"
        preds = pd.merge(preds, preds_m1, how="left", 
                         left_index=True, right_index=True)
        
        # Restrict to only certain atom types
        atom_set = {"H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"}
        preds = preds[["Res_name","Res_N","Res_type","Res_typem1"]+
                      list(atom_set.intersection(preds.columns))]
        
        preds.index = preds["Res_name"]
        
        self.preds = preds
        return(self.preds)
    
    
    def add_dummy_rows(self):
        """Add dummy rows to obs and preds to bring them to the same length.
        
        Also discard any atom types that aren't present in both obs and preds.
        """
        
        obs = self.obs.copy()
        preds = self.preds.copy()
        
        # Delete any prolines in preds
        preds = preds.drop(preds.index[preds["Res_type"]=="P"])
        
        # Restrict atom types
        # self.pars["atom_set"] is the set of atoms to be used in the analysis
        obs_metadata = list(set(obs.columns).difference(self.pars["atom_set"]))     
        preds_metadata = list(set(preds.columns).
                              difference(self.pars["atom_set"]))
        shared_atoms = list(self.pars["atom_set"].intersection(obs.columns).
                            intersection(preds.columns))
        obs = obs.loc[:,obs_metadata+shared_atoms]
        preds = preds.loc[:,preds_metadata+shared_atoms]
        
        # Create columns to keep track of dummy status
        preds["Dummy_res"] = False
        obs["Dummy_SS"] = False

        N = len(obs.index)
        M = len(preds.index)
        
        if N>M:     # If there are more spin systems than predictions
            dummies = pd.DataFrame(np.NaN, columns = preds.columns, 
                        index=["dummy_res_"+str(i) for i in 1+np.arange(N-M)])
            dummies["Res_name"] = dummies.index
            dummies["Dummy_res"] = True
            preds = preds.append(dummies)        
        elif M>N:
            dummies = pd.DataFrame(np.NaN, columns = obs.columns, 
                        index=["dummy_SS_"+str(i) for i in 1+np.arange(M-N)])
            dummies["SS_name"] = dummies.index
            dummies["Dummy_SS"] = True
            obs = obs.append(dummies)
            #obs.loc[["dummy_"+str(i) for i in 1+np.arange(M-N)]] = np.NaN
            #obs.loc[obs.index[N:M], "SS_name"] = ["dummy_"+str(i) for i in 1+np.arange(M-N)]

        self.obs = obs.copy()
        self.preds = preds.copy()
        
        return(self.obs, self.preds)
    
    def calc_log_prob_matrix(self, atom_sd=None, sf=1, default_prob=0.01, 
                             use_hadamac=False, cdf=False, rescale_delta=False, 
                             delta_correlation=False, shift_correlation=False,
                             verbose=False):
        """Calculate a matrix of -log10(match probabilities)
        
        use_hadamac: if True, amino acid type information will contribute to 
            the log probability
        cdf: if True, use cdf in probability matrix. Otherwise use pdf (cdf 
            uses chance of seeing a delta 'at least this great')
        rescale_delta: if True, the shift differences between obs and pred are 
            scaled so they are closer to the normal distribution
        delta_correlation: if True, correlated errors between different atom 
            types are accounted for in the probability 
        shift_correlation: if True, the correlation between observed shift and
            prediction error is accounted for.
        """
        
        # Use default atom_sd values if not defined
        if atom_sd==None:
            atom_sd = self.pars["atom_sd"]
#            atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231,
#                     'C':0.5330, 'CA':0.4412, 'CB':0.5163,
#                     'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}
        
        def calc_match_probability(obs, pred1):
            """ Calculate match scores between all observed spin systems and a 
            single predicted residue
            
            default_prob: probability assigned when an observation or 
                prediction is missing
            atom_sd: expected standard deviation for each atom type
            sf: scaling factor for the entire atom_sd dictionary
            use_hadamac: determines whether residue type information is used
            """
            
            # Throw away any non-atom columns
            obs_reduced = obs.loc[:, self.pars["atom_set"].
                                  intersection(obs.columns)]
            pred1_reduced = pred1.loc[self.pars["atom_set"].
                                      intersection(pred1.index)]
            
            # Calculate shift differences for each observed spin system
            delta = obs_reduced - pred1_reduced
            
            # Make a note of NA positions in delta, and set them to zero 
            # (this avoids warnings when using norm.cdf later)
            na_mask = delta.isna()
            delta[na_mask] = 0
            
            if self.pars["prob_method"] == "delta_correlation":
                overall_prob = pd.Series(index=delta.index)
                overall_prob[:] = 1
                
                d_mean = pd.read_csv("../data/d_mean.csv", header=None, 
                                     index_col=0).loc[delta.columns,1]
                d_cov = (pd.read_csv("../data/d_cov.csv", index_col=0).
                         loc[delta.columns,delta.columns])
                
                mvn = multivariate_normal(d_mean, d_cov)
                
                overall_prob = mvn.logpdf(delta)
                
                # Penalise missing shifts, unless also missing in predictions
                overall_prob = (overall_prob + log10(default_prob) * 
                            (na_mask.sum(axis=1) - pred1_reduced.isna().sum()))
                    
            else:
                prob = delta.copy()
                prob.iloc[:,:] = 1
                
                for c in delta.columns:
                    if self.pars["prob_method"] == "cdf":
                        # Use the cdf to calculate the probability of a 
                        # delta *at least* as great as the actual one
                        prob[c] = log10(2) + norm.logcdf(-1*abs(
                                pd.to_numeric(delta[c])), scale=atom_sd[c]*sf)
                    elif self.pars["prob_method"] == "pdf":
                        prob[c] = norm.logpdf(pd.to_numeric(delta[c]), 
                            scale=atom_sd[c]*sf)       
                    elif shift_correlation:
                        print("shift_correlation not yet implemented. Defaulting to pdf.")
                        prob[c] = norm.logpdf(pd.to_numeric(delta[c]), 
                            scale=atom_sd[c]*sf)
                    else:
                        print("Method for calculating probability not recognised. Defaulting to pdf.")
                        prob[c] = norm.logpdf(pd.to_numeric(delta[c]), 
                            scale=atom_sd[c]*sf)
                
                # In positions where data was missing, use default probability
                prob[na_mask] = log10(default_prob)
                
                # Calculate penalty for a HADAMAC mismatch
                if use_hadamac:
                    # If the i-1 aa type of the predicted residue matches the 
                    # HADAMAC group of the observation, probability is 1.
                    # Otherwise, probability defaults to 0.01
                    prob["HADAMACm1"] = 0.01
                    if type(pred1["Res_typem1"])==str:    # dummies have NaN
                        prob.loc[obs["HADAMACm1"].str.find(
                                pred1["Res_typem1"])>=0, "HADAMACm1"] = 1
            
                # Calculate overall probability of each row
                overall_prob = prob.sum(skipna=False, axis=1)
                
            return(overall_prob)
        
        obs = self.obs
        preds = self.preds
        
        # Initialise matrix as NaN
        log_prob_matrix = pd.DataFrame(np.NaN, index=obs.index, 
                                       columns=preds.index)    
        
        for i in preds.index:
            if verbose: print(i)
            log_prob_matrix.loc[:, i] = calc_match_probability(obs, 
                                                               preds.loc[i,:])
        
        
        # Calculate log of matrix
        log_prob_matrix[log_prob_matrix.isna()] = 2*np.nanmin(
                                                        log_prob_matrix.values)
        log_prob_matrix.loc[obs["Dummy_SS"], :] = 0
        log_prob_matrix.loc[:, preds["Dummy_res"]] = 0
        
        self.log_prob_matrix = log_prob_matrix
        return(self.log_prob_matrix)
        
    def calc_log_prob_matrix2(self, atom_sd=None, sf=1, default_prob=0.01, 
                             use_hadamac=False, cdf=False, 
                             delta_correlation=False, shift_correlation=False,
                             verbose=False):
        """Calculate a matrix of -log10(match probabilities)
        
        use_hadamac: if True, amino acid type information will contribute to 
            the log probability
        cdf: if True, use cdf in probability matrix. Otherwise use pdf (cdf 
            uses chance of seeing a delta 'at least this great')
        delta_correlation: if True, correlated errors between different atom 
            types are accounted for in the probability 
        shift_correlation: if True, the correlation between observed shift and
            prediction error is accounted for.
        """
        
        # Use default atom_sd values if not defined
        if atom_sd==None:
            atom_sd = self.pars["atom_sd"]
#            atom_sd={'H':0.1711, 'N':1.1169, 'HA':0.1231,
#                     'C':0.5330, 'CA':0.4412, 'CB':0.5163,
#                     'Cm1':0.5530, 'CAm1':0.4412, 'CBm1':0.5163}
        
        obs = self.obs
        preds = self.preds
        atoms = self.pars["atom_set"].intersection(obs.columns)
    
        log_prob_matrix = pd.DataFrame(0, index=obs.index, columns=preds.index)
        
        for atom in atoms:
            obs_atom = (obs[atom].repeat(len(obs.index)).values.
                        reshape([len(obs.index),-1]))
            preds_atom = (preds[atom].repeat(len(preds.index)).values.
                          reshape([len(preds.index),-1]).transpose())
            
            delta_atom = preds_atom - obs_atom
            
            # Make a note of NA positions in delta, and set them to zero 
            # (this avoids warnings when using norm.cdf later)
            na_mask = np.isnan(delta_atom)
            delta_atom[na_mask] = 0
            
            if self.pars["prob_method"] == "cdf":
                # Use the cdf to calculate the probability of a 
                # delta *at least* as great as the actual one
                prob_atom = pd.DataFrame(-2*norm.logcdf(abs(delta_atom), scale=atom_sd[atom]),
                                     index=obs.index, columns=preds.index)
            elif self.pars["prob_method"] == "pdf":
                prob_atom = pd.DataFrame(norm.logpdf(delta_atom, scale=atom_sd[atom]),
                                     index=obs.index, columns=preds.index)
            else:
                print("Method for calculating probability not recognised. Defaulting to pdf.")
                prob_atom = pd.DataFrame(norm.logpdf(delta_atom, scale=atom_sd[atom]),
                                     index=obs.index, columns=preds.index)
            
            
            prob_atom[na_mask] = log10(default_prob)
            
            log_prob_matrix = log_prob_matrix + prob_atom
    
        
        log_prob_matrix[log_prob_matrix.isna()] = 2*np.nanmin(
                                                        log_prob_matrix.values)
        log_prob_matrix.loc[obs["Dummy_SS"], :] = 0
        log_prob_matrix.loc[:, preds["Dummy_res"]] = 0
        
        self.log_prob_matrix = log_prob_matrix
        return(self.log_prob_matrix)
    
    def calc_dist_matrix(self, use_atoms=None, atom_scale=None, na_dist=0, rank=False):
        """Calculate the Euclidian distance between each observation and 
        prediction.
        
        use_atoms: limit the set of atoms considered. If None, uses 
            pars["atom_set"]
        atom_scale: how much the shift difference is scaled by
        na_dist: shift differences which can't be calculated (eg due to 
            missing data) are replaced with this value
        rank: if True, returns the rank of the distance per observation
        """
        obs = self.obs
        preds = self.preds
        
        # Use default atom_sd values if not defined
        if atom_scale==None:
            atom_scale = self.pars["atom_sd"]
        
        if use_atoms==None:
            atoms = self.pars["atom_set"].intersection(obs.columns)
        else:
            atoms = set(use_atoms).intersection(obs.columns)
        
        delta2 = pd.DataFrame(0, index=obs.index, columns=preds.index)
        
        for atom in atoms:
            obs_atom = (obs[atom].repeat(len(obs.index)).values.
                        reshape([len(obs.index),-1]))
            preds_atom = (preds[atom].repeat(len(preds.index)).values.
                          reshape([len(preds.index),-1]).transpose())
            
            delta2_atom = ((preds_atom - obs_atom)/atom_scale[atom])**2
            
            # Make a note of NA positions in delta, and set them to default value 
            na_mask = np.isnan(delta2_atom)
            delta2_atom[na_mask] = na_dist
            
            delta2 = delta2 + delta2_atom
            
        dist_matrix = delta2.applymap(sqrt)
        
        if rank:
            return (dist_matrix.rank(axis=1))
        else:
            return(dist_matrix)
    
    def find_best_assignment(self):
        """ Use the Hungarian algorithm to find the highest probability matching 
        (ie. the one with the lowest log probability sum)
        
        Return a data frame with matched observed and predicted shifts, and 
        the raw matching
        """
        
        obs = self.obs
        preds = self.preds
        log_prob_matrix = self.log_prob_matrix
        
        valid_atoms = list(self.pars["atom_set"])
        
        row_ind, col_ind = linear_sum_assignment(-1*log_prob_matrix)    
        # -1 because the algorithm minimises sum, but we want to maximise it.
        
        obs.index.name = "SS_index"
        preds.index.name = "Res_index"
        
        #Create assignment dataframe
        obs_names = [log_prob_matrix.index[r] for r in row_ind]
        pred_names = [log_prob_matrix.columns[c] for c in col_ind]
        
        assign_df = pd.DataFrame({
                "SS_name":obs_names,
                "Res_name":pred_names
                })
        
        # Merge residue info, shifts and predicteds into assignment dataframe
        assign_df = pd.merge(assign_df, 
                             preds[["Res_N","Res_type", "Res_name", 
                                    "Dummy_res"]], 
                             on="Res_name", how="outer")
        assign_df = assign_df[["Res_name","Res_N","Res_type","SS_name", 
                               "Dummy_res"]]
        assign_df = pd.merge(assign_df, 
                             obs.loc[:, obs.columns.isin(
                                     valid_atoms+["SS_name","Dummy_SS"])], 
                             on="SS_name", how="outer")
            # Above line raises an error about index/column confusion, 
            # which needs fixing. Now fixed?
        assign_df = pd.merge(assign_df, 
                             preds.loc[:, preds.columns.isin(
                                     valid_atoms+["Res_name"])],
                             on="Res_name", suffixes=("","_pred"), how="outer")
        
        assign_df["Log_prob"] = log_prob_matrix.lookup(
                                            log_prob_matrix.index[row_ind],
                                            log_prob_matrix.columns[col_ind])
        # Careful above not to get rows/columns confused
        assign_df = assign_df.sort_values(by="Res_N")
        
        self.assign_df = assign_df
        self.best_match_indexes = [row_ind, col_ind]
        return(assign_df, [row_ind, col_ind])
    
    def check_assignment_consistency(self, threshold=0.1):
        """ Find maximum mismatch and number of 'significant' mismatches for 
        each residue
        
        threshold: Minimum carbon shift difference for sequential residues to
            count as mismatched
        """
        
        assign_df = self.assign_df
        
        # First check if there are any sequential atoms
        carbons = pd.Series(["C","CA","CB"])
        carbons_m1 = carbons + "m1"
        seq_atoms = carbons[carbons.isin(assign_df.columns) & 
                            carbons_m1.isin(assign_df.columns)]
        seq_atoms_m1 = seq_atoms+"m1"
        #seq_atoms = list(seq_atoms)
    
        if seq_atoms.size==0:
            # You can't do a comparison
            assign_df["Max_mismatch_prev"] = np.NaN
            assign_df["Max_mismatch_next"] = np.NaN
            assign_df["Num_good_links_prev"] = np.NaN
            assign_df["Num_good_links_next"] = np.NaN
            return(assign_df)
        else:
            # First, get the i and i-1 shifts for the preceeding and 
            # succeeding residues
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
            tmp["Max_mismatch_prev"] = tmp["d"+seq_atoms+"_prev"].max(axis=1, 
                                                                   skipna=True)
            tmp["Max_mismatch_next"] = tmp["d"+seq_atoms+"_next"].max(axis=1,
                                                                   skipna=True)
            
            # Calculate number of consistent matches
            tmp["Num_good_links_prev"] = (tmp["d"+seq_atoms+"_prev"]<threshold).sum(axis=1)
            tmp["Num_good_links_next"] = (tmp["d"+seq_atoms+"_next"]<threshold).sum(axis=1)
            
            # Join relevant columns back onto assign_df
            tmp["Res_N"] = tmp.index
            assign_df = assign_df.join(tmp.loc[:,["Max_mismatch_prev", 
                                                  "Max_mismatch_next", 
                                                  "Num_good_links_prev", 
                                                  "Num_good_links_next"]], 
                                       on="Res_N")
            
            self.assign_df = assign_df
            return(self.assign_df)
    
    def find_alt_assignments(self, N=1, by_res=True,  verbose=False, 
                             return_full_assignments=False):
        """ Find the next-best assignment(s) for each residue or spin system
        
        This works by setting the log probability to a very high value for each 
        residue in turn, and rerunning the assignment
        
        Arguments:
        best_match_indexes: [row_ind, col_ind] output from find_best_assignment()
        N: number of alternative assignments to generate
        by_res: if true, calculate next best assignment for each residue. 
            Otherwise, calculate it for each spin system.
        
        Output:
        If return_full_assignments is False, returns a DataFrame containing the
        next-best assignment for each residue or spin system, and the relative
        overall probability of that assignment.
        If return_full_assignments is True, the output is a list with two 
        elements:
        1) the DataFrame of assignments and probabilities, as before
        2) a list of DataFrames containing the complete alternative assignments 
        for each residue.
        """
        
        obs = self.obs
        preds = self.preds
        log_prob_matrix = self.log_prob_matrix
        best_match_indexes = self.best_match_indexes
        
        # Calculate sum probability for the best matching
        best_sum_prob = sum(log_prob_matrix.lookup(
                              log_prob_matrix.index[best_match_indexes[0]],
                              log_prob_matrix.columns[best_match_indexes[1]]))
        
        penalty = 2*log_prob_matrix.min().min()     
        # This value is used to penalise the previous best match
        
        if by_res:
            # Create data structures for storing results
            assignment_list = [pd.Series(index = preds.index) for i in range(N)]
            rel_log_prob_list = [pd.Series(index = preds.index) for i in range(N)]
            results_dict = {}
            if return_full_assignments: 
                matching_df_list = [pd.DataFrame(index=preds.index, 
                                                 columns=preds.index) for i in range(N)]
            
            # Convert best_match_indexes to get a series of spin systems 
            # indexed by residue  
            best_matching = pd.Series(obs.index[best_match_indexes[0]], 
                                      index=preds.index[best_match_indexes[1]])
            alt_matching = pd.Series()
            
            for res in preds["Res_name"]:   # Consider each residue in turn
                if verbose: print(res)
                for i in range(N):
                    if i==0:
                        # Find the spin system that was matched to current 
                        # residue in optimal matching
                        ss = best_matching[res]
                        # Make a copy of the log_prob_matrix that can be safely 
                        # modified
                        tmp = log_prob_matrix.copy()    
                    else:
                        # Find the spin system that was matched to current
                        # residue in last round
                        ss = alt_matching[res]          
                    
                    # Penalise the match found for this residue
                    if obs.loc[ss, "Dummy_SS"]:
                        # If last match was a dummy spin system, penalise all 
                        # dummies
                        tmp.loc[obs["Dummy_SS"], res] = penalty     
                    else:
                        tmp.loc[ss, res] = penalty
                    
                    row_ind, col_ind = linear_sum_assignment(tmp)
                    
                    # Extract the info we want from the optimisation results
                    alt_matching = pd.Series(obs.index[row_ind], 
                                             index=preds.index[col_ind])
                    assignment_list[i][res] = alt_matching[res]
                    if return_full_assignments: 
                        matching_df_list[i].loc[:,res] = alt_matching
                    # Calculate the relative overall probability of this 
                    # assignment, compared to the optimal assignment
                    rel_log_prob_list[i][res] = sum(tmp.lookup(
                            alt_matching.values, 
                            alt_matching.index)) - best_sum_prob
        else:
            # Create data structures for storing results
            assignment_list = [pd.Series(index = obs.index) for i in range(N)]
            rel_log_prob_list = [pd.Series(index = obs.index) for i in range(N)]
            results_dict = {}
            if return_full_assignments: 
                matching_df_list = [pd.DataFrame(index=obs.index, columns=obs.index) for i in range(N)]
            
            # Convert best_match_indexes to get a series of spin systems 
            # indexed by spin system  
            best_matching = pd.Series(preds.index[best_match_indexes[1]], 
                                      index=obs.index[best_match_indexes[0]])
            alt_matching = pd.Series()
            
            for ss in obs["SS_name"]:   # Consider each spin system in turn
                if verbose: print(ss)
                for i in range(N):
                    if i==0:
                        # Find the residue that was matched to current spin 
                        # system in optimal matching
                        res = best_matching[ss]         
                        # Make a copy of the log_prob_matrix that can be safely 
                        # modified
                        tmp = log_prob_matrix.copy()    
                    else:
                        # Find the residue that was matched to current spin 
                        # system in last round
                        res = alt_matching[ss]          
                    
                    # Penalise the match found for this residue
                    if preds.loc[res, "Dummy_res"]:
                        # If last match was a dummy residue, penalise all 
                        # dummies
                        tmp.loc[ss, preds["Dummy_res"]] = penalty     
                    else:
                        tmp.loc[ss, res] = penalty
                    
                    row_ind, col_ind = linear_sum_assignment(tmp)
                    
                    # Extract the info we want from the optimisation results
                    alt_matching = pd.Series(preds.index[col_ind], 
                                             index=obs.index[row_ind])
                    assignment_list[i][ss] = alt_matching[ss]
                    if return_full_assignments: 
                        matching_df_list[i].loc[:,ss] = alt_matching
                    # Calculate the relative overall probability of this 
                    # assignment, compared to the optimal assignment
                    rel_log_prob_list[i][ss] = sum(tmp.lookup(
                            alt_matching.index, 
                            alt_matching.values)) - best_sum_prob            
            
        # Store the results as a dataframe
        for i in range(N):
            results_dict["Alt_assign_"+str(i+1)] = assignment_list[i]
            results_dict["Alt_rel_prob"+str(i+1)] = rel_log_prob_list[i]
        results_df = pd.DataFrame(results_dict)
        if by_res:
            results_df["Res_name"] = results_df.index
        else:
            results_df["SS_name"] = results_df.index
            
        if return_full_assignments:
            return(results_df, matching_df_list)
        else:
            return(results_df)
            
    def find_alt_assignments2(self, N=1, by_ss=True, verbose=False):
        """ Find the next-best assignment(s) for each residue or spin system
        
        This works by setting the log probability to a very high value for each 
        residue in turn, and rerunning the assignment
        
        Arguments:
        best_match_indexes: [row_ind, col_ind] output from find_best_assignment()
        N: number of alternative assignments to generate
        by_ss: if true, calculate next best assignment for each spin system. 
            Otherwise, calculate it for each residue.
        
        Output:
        A Dataframe containing the original assignments, and the 
        alt_assignments by
        """
        
        obs = self.obs
        preds = self.preds
        log_prob_matrix = self.log_prob_matrix
        best_match_indexes = self.best_match_indexes
        
        # Calculate sum probability for the best matching
        best_sum_prob = sum(log_prob_matrix.lookup(
                log_prob_matrix.index[best_match_indexes[0]],
                log_prob_matrix.columns[best_match_indexes[1]]))
        
        # Calculate the value used to penalise the best match for each residue
        penalty = 2*log_prob_matrix.min().min()     
        logging.debug("Penalty value: %f", penalty)
        
        # Initialise DataFrame for storing alt_assignments
        self.alt_assign_df = self.assign_df.copy()
        self.alt_assign_df["Rank"] = 1
        self.alt_assign_df["Rel_prob"] = 0
       
        if by_ss:
            # Convert best_match_indexes to get a series of spin systems 
            # indexed by spin system  
            best_matching = pd.Series(preds.index[best_match_indexes[1]], 
                                      index=obs.index[best_match_indexes[0]])
            
            for ss in obs["SS_name"]:   # Consider each spin system in turn
                logging.debug("Finding alt assignments for spin system %s", ss)
                if verbose: print(ss)
                a = deepcopy(self)
                for i in range(N):
                    if i==0:
                        # Find the residue that was matched to current spin 
                        # system in optimal matching
                        res = best_matching[ss]         
                    else:
                        # Find the residue that was matched to current spin 
                        # system in last round
                        res = alt_match          
                    
                    # Penalise the match found for this residue
                    if preds.loc[res, "Dummy_res"]:
                        a.log_prob_matrix.loc[ss, preds["Dummy_res"]] = penalty     # If last match was a dummy residue, penalise all dummies
                    else:
                        a.log_prob_matrix.loc[ss, res] = penalty
                    
                    a.find_best_assignment()
                    a.assign_df.index = a.assign_df["SS_name"]
                    a.assign_df["Rank"] = i+2
                    alt_sum_prob = sum(a.log_prob_matrix.lookup(
                            a.log_prob_matrix.index[a.best_match_indexes[0]],
                            a.log_prob_matrix.columns[a.best_match_indexes[1]]))
                    a.assign_df["Rel_prob"] = alt_sum_prob - best_sum_prob
                    
                    alt_match = a.assign_df.loc[ss, "Res_name"] 
                    self.alt_assign_df = self.alt_assign_df.append(
                                    a.assign_df.loc[ss, :], ignore_index=True)
                self.alt_assign_df = self.alt_assign_df.sort_values(
                                                        by=["SS_name", "Rank"])
                
        else:
            # Convert best_match_indexes to get a series of spin systems 
            # indexed by spin system  
            best_matching = pd.Series(obs.index[best_match_indexes[0]], 
                                      index=preds.index[best_match_indexes[1]])
            
            for res in preds["Res_name"]:   # Consider each spin system in turn
                if verbose: print(res)
                a = deepcopy(self)
                for i in range(N):
                    if i==0:
                        # Find the spin system that was matched to current 
                        # residue in optimal matching
                        ss = best_matching[res]         
                    else:
                        # Find the spin system that was matched to current 
                        # residue in last round
                        ss = alt_match          
                    
                    # Penalise the match found for this residue
                    if obs.loc[ss, "Dummy_SS"]:
                        # If last match was a dummy residue, penalise all dummies
                        a.log_prob_matrix.loc[obs["Dummy_SS"], res] = penalty     
                    else:
                        a.log_prob_matrix.loc[ss, res] = penalty
                    
                    a.find_best_assignment()
                    a.assign_df.index = a.assign_df["Res_name"]
                    a.assign_df["Rank"] = i+2
                    alt_sum_prob = sum(a.log_prob_matrix.lookup(
                            a.log_prob_matrix.index[a.best_match_indexes[0]],
                            a.log_prob_matrix.columns[a.best_match_indexes[1]]))
                    a.assign_df["Rel_prob"] = alt_sum_prob - best_sum_prob
                    
                    alt_match = a.assign_df.loc[res, "SS_name"] 
                    self.alt_assign_df = self.alt_assign_df.append(
                                    a.assign_df.loc[res, :], ignore_index=True)
                self.alt_assign_df = self.alt_assign_df.sort_values(
                                                    by=["Res_name", "Rank"])
        
        return(self.alt_assign_df)

    def output_peaklists(self, filepath, format="sparky", 
                         spectra=["hsqc","hnco","hncaco","hncacb", "hncocacb"]):
        """ Output assigned peaklists
        """
        return(0)
    
    def plot_strips(self, atom_list=["C","Cm1","CA","CAm1","CB","CBm1"]):
        """ Make a strip plot of the assignment
        
        atom_list: only plot data for these atom types
        """
        assign_df = self.assign_df
        
        # Narrow down atom list to those actually present
        atom_list = list(set(atom_list).intersection(assign_df.columns))
        
        # First, convert assign_df from wide to long
        plot_df = assign_df.loc[:,["Res_N", "Res_type", "Res_name", "SS_name", 
                                   "Dummy_res", "Dummy_SS"]+atom_list]
        plot_df = plot_df.melt(id_vars=["Res_N", "Res_type", "Res_name", 
                                        "SS_name", "Dummy_res", "Dummy_SS"],
                                   value_vars=atom_list, var_name="Atom_type",
                                   value_name="Shift")
        
        # Add columns with information to be plotted
        plot_df["i"] = "0"     # Track if shift is from the i or i-1 residue
        plot_df.loc[plot_df["Atom_type"].isin(["Cm1","CAm1","CBm1"]),"i"] = "-1"
        plot_df["Atom_type"] = plot_df["Atom_type"].replace({"Cm1":"C", 
                                                   "CAm1":"CA", "CBm1":"CB"}) 
                                                    # Simplify atom type
        
        plot_df["seq_group"] = plot_df["Res_N"] + plot_df["i"].astype("int")
        
        # Pad Res_name column with spaces so that sorting works correctly
        plot_df["Res_name"] = plot_df["Res_name"].str.pad(6)
        plot_df["x_name"] = plot_df["Res_name"] + "_(" + plot_df["SS_name"] + ")"
        
        # Make the plot
        plt = ggplot(aes(x="x_name"), data=plot_df) 
        plt = plt + geom_point(aes(y="Shift", colour="i", shape="Dummy_res"))
        plt = plt + scale_y_reverse() + scale_shape_manual(values=["o","x"])
        # Add lines connecting i to i-1 points
        plt = plt + geom_line(aes(y="Shift", group="seq_group"), 
                              data=plot_df.loc[~plot_df["Dummy_res"],])        
        plt = plt + geom_line(aes(y="Shift", group="x_name"), linetype="dashed")
        plt = plt + facet_grid("Atom_type~.", scales="free") 
        plt = plt + scale_colour_brewer(type="Qualitative", palette="Set1") 
        plt = plt + xlab("Residue name") + ylab("Chemical shift (ppm)")
        plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90))
        
        return(plt)
    
    def plot_seq_mismatch(self):
        """ Make a plot of the maximum sequential mismatch between i-1, i and 
        i+1 residues
        """
        assign_df = self.assign_df
        
        # Check that the assignment data frame has the right columns
        if not all(pd.Series(['Max_mismatch_prev', 'Max_mismatch_next']).
                   isin(assign_df.columns)):
            return(None)
        else:
            # Pad Res_name column with spaces so that sorting works correctly
            assign_df["Res_name"] = assign_df["Res_name"].str.pad(6)
            assign_df["x_name"] = (assign_df["Res_name"] + "_(" + 
                                     assign_df["SS_name"] + ")")
            
            # Make the plot
            plt = ggplot(aes(x="x_name"), data=assign_df) 
            plt = plt + geom_col(aes(y="abs(Max_mismatch_prev)"))
            plt = plt + xlab("Residue name")
            plt = plt + ylab("Mismatch to previous residue (ppm)")
            plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90))
                   
            return(plt)

#%%
        
#### Testing 

#a = NAPS_assigner()
#a.import_obs_shifts("~/GitHub/NAPS/data/testset/simplified_BMRB/6338.txt")
#a.import_pred_shifts("~/GitHub/NAPS/data/P3a_L273R/shiftx2.cs", offset=208)
#a.add_dummy_rows()
#a.calc_log_prob_matrix(sf=2, verbose=False)
#assign_df, best_match_indexes = a.find_best_assignment()


## Import the observed and predicted shifts
#    obs = import_obs_shifts(obs_file)
#    preds = import_pred_shifts(preds_file)
#    
#    # Add dummy rows so that obs and preds are the same length
#    obs, preds = add_dummy_rows(obs, preds)
#    
#    # Calculate the log probability for each observation/prediction pair
#    log_prob_matrix = calc_log_prob_matrix(obs, preds, sf=2)
#    
#    # Find the assignment with the highest overall probability
#    assign_df, matching = find_best_assignment(obs, preds, log_prob_matrix)        
        
        