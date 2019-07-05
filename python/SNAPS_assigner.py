#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Defines a class with functions related to the actual assignment process, and 
outputting the results.

@author: aph516
"""

import numpy as np
import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
from bokeh.models.ranges import Range1d
from bokeh.models import WheelZoomTool, LabelSet, ColumnDataSource, Span
from bokeh.io import export_png
from bokeh.embed import json_item
from scipy.stats import norm, multivariate_normal
from scipy.optimize import linear_sum_assignment
from math import log10, sqrt
from copy import deepcopy
from pathlib import Path
#from Bio.SeqUtils import seq1
from distutils.util import strtobool
from collections import namedtuple
from sortedcontainers import SortedListWithKey
import logging
import yaml

class SNAPS_assigner:
    # Functions
    def __init__(self):
        self.obs = None
        self.preds = None
        self.log_prob_matrix = None
        self.mismatch_matrix = None
        self.consistent_links_matrix = None
        self.assign_df = None
        self.alt_assign_df = None
        self.best_match_indexes = None
        self.pars = {"pred_offset": 0,
                "prob_method": "pdf",
                "pred_correction": False,
                "delta_correlation": False,
                "alt_assignments": 1,
                "atom_set": {"H","N","HA","C","CA","CB","C_m1","CA_m1","CB_m1"},
                "atom_sd": {'H':0.1711, 'N':1.1169, 'HA':0.1231,
                            'C':0.5330, 'CA':0.4412, 'CB':0.5163,
                            'C_m1':0.5530, 'CA_m1':0.4412, 'CB_m1':0.5163},
                "plot_strips": False,
                "seq_link_threshold": 0.1}
            
    def read_config_file(self, filename):
        "Read a configuration file written in YAML format"       
        f = open(filename, 'r')
        
        self.pars = yaml.safe_load(f)
        self.pars["atom_set"] = set(self.pars["atom_set"])
        
        # TODO: Add a check that all necessary parameters are present
        
        return(self.pars)
    
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
                        break
                        
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
            
        # Convert from long to wide format
        preds = preds_long.pivot(index="Res_N", columns="Atom_type", 
                                 values="Shift")
        
        # Add the other data back in
        tmp = preds_long[["Res_N","Res_type","Res_name"]]
        tmp = tmp.drop_duplicates(subset="Res_name")
        tmp.index = tmp["Res_N"]
        tmp.index.name = None
        preds = pd.concat([tmp, preds], axis=1)
        
        # Make columns for the i-1 predicted shifts of C, CA and CB
        preds_m1 = preds[list({"C","CA","CB","Res_type","Res_name"}.
                              intersection(preds.columns))].copy()
        preds_m1.index = preds_m1.index+1
        preds_m1.columns = preds_m1.columns + "_m1"
        preds = pd.merge(preds, preds_m1, how="left", 
                         left_index=True, right_index=True)
        
        # Make column for the i+1 Res_name
        preds_p1 = preds[["Res_name"]].copy()
        preds_p1.index = preds_p1.index-1
        preds_p1.columns = ["Res_name_p1"]
        preds = pd.merge(preds, preds_p1, how="left", 
                         left_index=True, right_index=True)
        
        # Restrict to only certain atom types
        atom_set = {"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}
        preds = preds[["Res_name","Res_N","Res_type","Res_name_m1","Res_name_p1","Res_type_m1"]+
                      list(atom_set.intersection(preds.columns))]
        
        self.preds = preds
        return(self.preds)
    
    def simulate_pred_shifts(self, filename, sd, seed=None):
        """Generate a 'simulated' predicted shift DataFrame by importing some 
        observed chemical shifts (in 'test' format), and adding Gaussian 
        errors.
        
        sd: dictionary of standard deviation for each atom type 
            eg. {"H":0.1,"N":0.5}
        """
        from SNAPS_importer import SNAPS_importer
        import numpy.random as rand
        
        if seed is not None:
            rand.seed(seed)
            
        
        importer = SNAPS_importer()
        importer.import_testset_shifts(filename, remove_Pro=False)
        
        preds = importer.obs
        
        # Limit columns
        preds = preds[["Res_N","Res_type"]+
                      list({"C","CA","CB","H","N","HA"}.intersection(preds.columns))]
        
        # Add the random shifts
        for atom in self.pars["atom_set"].intersection(preds.columns):
            preds.loc[:,atom] += rand.normal(0, sd[atom], size=len(preds.index))
        
        # Add other columns back in
        preds.insert(1, "Res_name", (preds["Res_N"].astype(str) + 
                  preds["Res_type"]))
        preds.loc[:,"Res_name"] = [s.rjust(5) for s in preds["Res_name"]]
        
        preds.index = preds["Res_N"]
        preds.index.name=None
        
        # Make columns for the i-1 predicted shifts of C, CA and CB
        preds_m1 = preds[list({"C","CA","CB","Res_type","Res_name"}.
                              intersection(preds.columns))].copy()
        preds_m1.index = preds_m1.index+1
        preds_m1.columns = preds_m1.columns + "_m1"
        preds = pd.merge(preds, preds_m1, how="left", 
                         left_index=True, right_index=True)
        
        # Make column for the i+1 Res_name
        preds_p1 = preds[["Res_name"]].copy()
        preds_p1.index = preds_p1.index-1
        preds_p1.columns = ["Res_name_p1"]
        preds = pd.merge(preds, preds_p1, how="left", 
                         left_index=True, right_index=True)
        
        # Restrict to only certain atom types
        atom_set = {"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}
        preds = preds[["Res_name","Res_N","Res_type","Res_name_m1","Res_name_p1","Res_type_m1"]+
                      list(atom_set.intersection(preds.columns))]
        
        self.preds = preds
        return(self.preds)
    
    def add_dummy_rows(self):
        """Add dummy rows to obs and preds to bring them to the same length.
        
        Also discard any atom types that aren't present in both obs and preds.
        """
        
        obs = self.obs.copy()
        preds = self.preds.copy()
        
        # Delete any prolines in preds
        self.all_preds = preds.copy()   # Keep a copy of all predictions
        preds = preds.drop(preds.index[preds["Res_type"]=="P"])
        
        # Create a data frame to quickly look up the i-1 and i+1 Res_name
        preds.index = preds["Res_name"]
        preds.index.name = None
        self.neighbour_df = preds[["Res_name_m1","Res_name_p1"]]
        #self.Res_name_m1.loc[preds["Res_type_m1"]=="P"] = np.NaN
        
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
                        index=["DR_"+str(i) for i in 1+np.arange(N-M)])
            dummies["Res_name"] = dummies.index
            dummies["Dummy_res"] = True
            preds = preds.append(dummies)        
        elif M>N:
            dummies = pd.DataFrame(np.NaN, columns = obs.columns, 
                        index=["DSS_"+str(i) for i in 1+np.arange(M-N)])
            dummies["SS_name"] = dummies.index
            dummies["Dummy_SS"] = True
            obs = obs.append(dummies)
            #obs.loc[["dummy_"+str(i) for i in 1+np.arange(M-N)]] = np.NaN
            #obs.loc[obs.index[N:M], "SS_name"] = ["dummy_"+str(i) for i in 1+np.arange(M-N)]
        
        # Set missing Res_name_m1 entries to NaN
        self.neighbour_df.loc[~self.neighbour_df["Res_name_m1"].isin(preds["Res_name"]), 
                          "Res_name_m1"] = np.NaN
        self.neighbour_df.loc[~self.neighbour_df["Res_name_p1"].isin(preds["Res_name"]), 
                          "Res_name_p1"] = np.NaN
        
        self.obs = obs.copy()
        self.preds = preds.copy()
        
        return(self.obs, self.preds)
    
        
    def calc_log_prob_matrix(self, atom_sd=None, sf=1, default_prob=0.01, 
                             verbose=False):
        """Calculate a matrix of -log10(match probabilities)
        
        sf: Multiply the provided atom_sd's by this number
        default_prob: penalty for missing data
        """
        
        # Use default atom_sd values if not defined
        if atom_sd==None:
            atom_sd = self.pars["atom_sd"]
        
        obs = self.obs
        preds = self.preds
        atoms = list(self.pars["atom_set"].intersection(obs.columns))
        
        if self.pars["pred_correction"]:
            # Import parameters for correcting the shifts
            lm_pars = pd.read_csv(self.pars["pred_correction_file"], index_col=0)
        
        if self.pars["delta_correlation"]:
            # Import parameters describing the delta correlations
            if self.pars["pred_correction"]:
                d_mean = pd.read_csv(self.pars["delta_correlation_mean_corrected_file"], 
                                     header=None, index_col=0).loc[atoms,1]
                d_cov = (pd.read_csv(self.pars["delta_correlation_cov_corrected_file"], 
                                     index_col=0).loc[atoms,atoms])
            else:
                d_mean = pd.read_csv(self.pars["delta_correlation_mean_file"], 
                                     header=None, index_col=0).loc[atoms,1]
                d_cov = (pd.read_csv(self.pars["delta_correlation_cov_file"], 
                                     index_col=0).loc[atoms,atoms])
            delta_list = []
        
                   
        
        log_prob_matrix = pd.DataFrame(0, index=obs.index, columns=preds.index)
        for atom in atoms:
            # The most efficient way I've found to do the calculation is to 
            # take the obs and preds shift columns for an atom, repeat each 
            # into a matrix, then subtract these matrixes from each other. 
            # That way, all calculations take advantage of vectorisation. 
            # Much faster than using loops.
            obs_atom = pd.DataFrame(obs[atom].repeat(len(obs.index)).values.
                                reshape([len(obs.index),-1]),
                                index=obs.index, columns=preds.index)
            preds_atom = pd.DataFrame(preds[atom].repeat(len(preds.index)).values.
                                reshape([len(preds.index),-1]).transpose(),
                                index=obs.index, columns=preds.index)
            
            
            # If predicting corrections, apply a linear transformation of delta
            if self.pars["pred_correction"]:
                preds_corr_atom = preds_atom
                for res in preds["Res_type"].dropna().unique():
                    if (atom+"_"+res) in lm_pars.index:
                         
                        grad = lm_pars.loc[(lm_pars["Atom_type"]==atom) & 
                                           (lm_pars["Res_type"]==res),
                                           "Grad"].tolist()[0]
                        offset = lm_pars.loc[(lm_pars["Atom_type"]==atom) & 
                                             (lm_pars["Res_type"]==res),
                                             "Offset"].tolist()[0]
                        if atom in ("C_m1","CA_m1","CB_m1"):
                            preds_corr_atom.loc[:,preds["Res_type_m1"]==res] = (
                                    preds_atom.loc[:, preds["Res_type_m1"]==res]
                                    - grad * obs_atom.loc[:, preds["Res_type_m1"]==res]
                                    - offset) 
                        else:
                            preds_corr_atom.loc[:, preds["Res_type"]==res] = (
                                    preds_atom.loc[:, preds["Res_type"]==res]
                                    - grad * obs_atom.loc[:, preds["Res_type"]==res]
                                    - offset)
                                
                delta_atom = preds_corr_atom - obs_atom
            else:
                delta_atom = preds_atom - obs_atom
            
            if self.pars["delta_correlation"]:
                delta_list = delta_list + [delta_atom.values]
            else:
                # Make a note of NA positions in delta, and set them to zero 
                # (this avoids warnings when using norm.cdf later)
                na_mask = np.isnan(delta_atom)
                delta_atom[na_mask] = 0
                
                if self.pars["prob_method"] == "cdf":
                    # Use the cdf to calculate the probability of a 
                    # delta *at least* as great as the actual one
                    prob_atom = pd.DataFrame(-2*norm.logcdf(abs(delta_atom), 
                                                            scale=atom_sd[atom]),
                                             index=obs.index, columns=preds.index)
                elif self.pars["prob_method"] == "pdf":
                    prob_atom = pd.DataFrame(norm.logpdf(delta_atom, 
                                                         scale=atom_sd[atom]),
                                             index=obs.index, columns=preds.index)
                else:
                    print("Method for calculating probability not recognised. Defaulting to pdf.")
                    prob_atom = pd.DataFrame(norm.logpdf(delta_atom, scale=atom_sd[atom]),
                                         index=obs.index, columns=preds.index)
                
                
                prob_atom[na_mask] = log10(default_prob)
                
                log_prob_matrix = log_prob_matrix + prob_atom
        
        if self.pars["delta_correlation"]:
            delta_mat = np.array(delta_list)
            delta_mat = np.moveaxis(delta_mat, 0, -1)
            #self.delta_mat = delta_mat
            
            # Make a note of NA positions in delta, and set them to zero
            na_mask = np.isnan(delta_mat)
            delta_mat[na_mask] = 0
            
            # TODO: Check that atom order in list is the same as in d_mean, d_cov
            mvn = multivariate_normal(d_mean, d_cov)
            log_prob_matrix = pd.DataFrame(mvn.logpdf(delta_mat), 
                                           index=obs.index, columns=preds.index)
            
            # TODO: Need to apply correction for missing data?
        
        if self.pars["use_ss_class_info"]:
            # For each type of residue type information that's available, make a 
            # matrix showing the probability modifications due to type mismatch, 
            # then add it to log_prob_matrix
            # Maybe make SS_class mismatch a parameter in config file?
            for ss_class in {"SS_class","SS_class_m1"}.intersection(obs.columns):
                #print(ss_class)
                SS_class_matrix = pd.DataFrame(0, index=log_prob_matrix.index, 
                                           columns=log_prob_matrix.columns)
                
                # For each amino acid type in turn:
                for res in preds["Res_type"].dropna().unique():                
                    # Work out which observations could be that aa type
                    allowed = obs[ss_class].str.contains(res).fillna(True)
                    # Select the predictions which are that aa type
                    pred_list = preds.loc[preds["Res_type_m1"]==res,"Res_name"]
                    # For the selected predictions, penalise any observations 
                    # where the current aa type is not allowed
                    for p in pred_list:
                        SS_class_matrix.loc[:,p] = (~allowed)*-100 #log10(0.01)
            
                log_prob_matrix = log_prob_matrix + SS_class_matrix
                
            
        
        log_prob_matrix[log_prob_matrix.isna()] = 2*np.nanmin(
                                                        log_prob_matrix.values)
        log_prob_matrix.loc[obs["Dummy_SS"], :] = 0
        log_prob_matrix.loc[:, preds["Dummy_res"]] = 0
        
        self.log_prob_matrix = log_prob_matrix
        return(self.log_prob_matrix)
    
    def find_best_assignments(self, inc=None, exc=None, return_none_if_all_dummy=False, verbose=True):
        """ Use the Hungarian algorithm to find the highest probability matching 
        (ie. the one with the lowest log probability sum), with constraints.
        
        Returns a data frame with the SS_names and Res_names of the matching. 
        (Doesn't change the internal state of the SNAPS_assigner instance.)
        
        inc: a DataFrame of (SS,Res) pairs which must be part of the assignment. 
            First column has the SS_names, second has the Res_names .
        exc: a DataFrame of (SS,Res) pairs which may not be part of the assignment.
        return_none_if_all_dummy: Sometimes after setting aside the pairs that 
            must be included in the final assignment, only dummy residues or 
            spin systems will remain. If this argument is True, the function 
            will return None in these cases.
        verbose: if True, prints messages if inc and exc contain reduntant info
        """
        obs = self.obs
        preds = self.preds
        log_prob_matrix = deepcopy(self.log_prob_matrix)
        
        if inc is not None:
            # Check for conflicting entries in inc
            conflicts = inc["SS_name"].duplicated(keep=False) | inc["Res_name"].duplicated(keep=False)
            if any(conflicts):
                print("Error: entries in inc conflict with one another - dropping conflicts")
                print(inc[conflicts])
                inc = inc[~conflicts]
            
            if exc is not None:
                # Check constraints are consistent
                # Get rid of any entries in exc which share a Res or SS with inc
                exc_in_inc = exc["SS_name"].isin(inc["SS_name"]) | exc["Res_name"].isin(inc["Res_name"])
                if any(exc_in_inc):
                    if verbose:
                        print("Some values in exc are also found in inc, so are redundant.")
                        print(exc[exc_in_inc])
                    exc = exc.loc[~exc_in_inc, :]
                    
            # Removed fixed assignments from probability matrix and obs, preds 
            # dataframes. Latter is needed if inc includes any dummy SS/res, 
            # and to detect if the reduced data is entirely dummies
            log_prob_matrix_reduced = log_prob_matrix.drop(index=inc["SS_name"]).drop(columns=inc["Res_name"])
            obs_reduced = obs.drop(index=inc["SS_name"])
            preds_reduced = preds.drop(index=inc["Res_name"])
        else:
            log_prob_matrix_reduced = log_prob_matrix
            obs_reduced = obs
            preds_reduced = preds
        
        if return_none_if_all_dummy:
            if obs_reduced["Dummy_SS"].all() or preds_reduced["Dummy_res"].all():
                return(None)
        
        if exc is not None:
            # Penalise excluded SS,Res pairs
            penalty = 2*log_prob_matrix.min().min()
            for index, row in exc.iterrows():
                # Need to account for dummy residues or spin systems
                if preds.loc[row["Res_name"], "Dummy_res"]:
                    log_prob_matrix_reduced.loc[row["SS_name"], 
                                                preds_reduced["Dummy_res"]] = penalty
                elif obs.loc[row["SS_name"], "Dummy_SS"]:
                    log_prob_matrix_reduced.loc[obs_reduced["Dummy_SS"], 
                                                row["Res_name"]] = penalty
                else:
                    log_prob_matrix_reduced.loc[row["SS_name"], row["Res_name"]] = penalty
        
        row_ind, col_ind = linear_sum_assignment(-1*log_prob_matrix_reduced)
        # -1 because the algorithm minimises sum, but we want to maximise it.
        
        # Construct results dataframe
        matching_reduced = pd.DataFrame({"SS_name":log_prob_matrix_reduced.index[row_ind],
                                           "Res_name":log_prob_matrix_reduced.columns[col_ind]})
        
        if inc is not None:
            matching = pd.concat([inc, matching_reduced])             
            return(matching)
        else:
            return(matching_reduced)
    
    def make_assign_df(self, matching, set_assign_df=False):
        """Make a dataframe with full assignment information, given a dataframe 
        of SS_name and Res_name.
        
        Matching may have additional columns, which will also be kept.
        """
        obs = self.obs
        preds = self.preds
        log_prob_matrix = self.log_prob_matrix
        valid_atoms = list(self.pars["atom_set"])
        extra_cols = set(matching.columns).difference({"SS_name","Res_name"})
        
        obs.index.name = None
        
        assign_df = pd.merge(matching, 
                             preds.loc[:,["Res_N","Res_type", "Res_name", 
                                    "Dummy_res"]], 
                             on="Res_name", how="left")
        assign_df = assign_df[["Res_name","Res_N","Res_type","SS_name", 
                               "Dummy_res"]+list(extra_cols)]
        assign_df = pd.merge(assign_df, 
                             obs.loc[:, obs.columns.isin(
                                     ["SS_name","Dummy_SS"]+valid_atoms)], 
                             on="SS_name", how="left")
        assign_df = pd.merge(assign_df, 
                             preds.loc[:, preds.columns.isin(
                                     valid_atoms+["Res_name"])],
                             on="Res_name", suffixes=("","_pred"), how="left")
        
        assign_df["Log_prob"] = log_prob_matrix.lookup(
                                            assign_df["SS_name"],
                                            assign_df["Res_name"])
        # Careful above not to get rows/columns confused
        
        assign_df = assign_df.sort_values(by="Res_N")
        
        if set_assign_df:
            self.assign_df = assign_df
        
        return(assign_df)
    
    def calc_overall_matching_prob(self, matching):
        "Calculate the sum log probability of a particular matching"
        return(sum(self.log_prob_matrix.lookup(matching["SS_name"], 
                                                   matching["Res_name"])))    
    
    def calc_mismatch_matrix(self, threshold=0.2):
        """Calculate matrix of the mismatch between i and i-1 observed shifts 
        for all spin system pairs. Also, make matrix of number of consistent 
        links for all spin system pairs.
        """
        obs = self.obs
        
        # First check if there are any sequential atoms
        carbons = pd.Series(["C","CA","CB"])
        carbons_m1 = carbons + "_m1"
        seq_atoms = carbons[carbons.isin(obs.columns) & 
                            carbons_m1.isin(obs.columns)]
        #seq_atoms_m1 = seq_atoms+"_m1"
        #seq_atoms = list(seq_atoms)
    
        if seq_atoms.size==0:
            # You can't make a mismatch matrix
            return(None)
        else:
            mismatch_matrix = pd.DataFrame(0, index=obs["SS_name"], columns=obs["SS_name"])
            mismatch_matrix.columns.name = "i"
            mismatch_matrix.index.name = "i_m1"
            
            consistent_links_matrix = mismatch_matrix.copy()
            
            for atom in seq_atoms:
                i_atom = (obs[atom].repeat(len(obs.index)).values.
                          reshape([len(obs.index),-1]))
                i_m1_atom = (obs[atom+"_m1"].repeat(len(obs.index)).values.
                          reshape([len(obs.index),-1]).transpose())
                mismatch_atom = pd.DataFrame(i_atom - i_m1_atom)
                mismatch_atom.columns = obs["SS_name"]
                mismatch_atom.columns.name = "i"
                mismatch_atom.index = obs["SS_name"]
                mismatch_atom.index.name = "i_m1"
                
                consistent_links_atom = (abs(mismatch_atom) < threshold)
                
                # Make a note of NA positions in delta, and set them to default value 
                na_mask = np.isnan(mismatch_atom)
                mismatch_atom[na_mask] = 0
                consistent_links_atom[na_mask] = 0
                
                # Update mismatch matrix
                mismatch_matrix = mismatch_matrix.combine(abs(mismatch_atom), np.maximum)
                consistent_links_matrix = consistent_links_matrix + consistent_links_atom
            
            self.mismatch_matrix = mismatch_matrix
            self.consistent_links_matrix = consistent_links_matrix
            return(self.mismatch_matrix)
    
    def check_matching_consistency(self, matching, threshold=0.2):
        """Calculate mismatch scores and assignment confidence for a given matching"""
        # Create mismatch matrix if it doesn't already exist
        if self.mismatch_matrix is None:
            self.calc_mismatch_matrix()
        
        
        # Add a Res_name_m1 column to matching DataFrame
        matching.index = matching["Res_name"]
        matching.index.name = None
        
        tmp = pd.concat([matching, self.neighbour_df], axis=1)
        
        # Add a SS_name_m1 and SS_name_p1 columns
        tmp = tmp.merge(matching[["SS_name"]], how="left", left_on="Res_name_m1", 
                        right_index=True, suffixes=("","_m1"))
        tmp = tmp.merge(matching[["SS_name"]], how="left", left_on="Res_name_p1", 
                        right_index=True, suffixes=("","_p1"))
        
        # Filter out rows with NaN's in SS_name_m1/p1
        tmp_m1 = tmp.dropna(subset=["SS_name_m1"])
        tmp_p1 = tmp.dropna(subset=["SS_name_p1"])
        
        # Get the mismatches
        tmp["Max_mismatch_m1"] = pd.Series(self.mismatch_matrix.lookup(
                                tmp_m1["SS_name_m1"], tmp_m1["SS_name"]),
                                index=tmp_m1.index)
        tmp["Max_mismatch_p1"] = pd.Series(self.mismatch_matrix.lookup(
                                tmp_p1["SS_name"], tmp_p1["SS_name_p1"]),
                                index=tmp_p1.index)
        tmp["Max_mismatch"] = tmp[["Max_mismatch_m1","Max_mismatch_p1"]].max(axis=1)
        
        tmp["Num_good_links_m1"] = pd.Series(self.consistent_links_matrix.lookup(
                                tmp_m1["SS_name_m1"], tmp_m1["SS_name"]),
                                index=tmp_m1.index)
        tmp["Num_good_links_p1"] = pd.Series(self.consistent_links_matrix.lookup(
                                tmp_p1["SS_name"], tmp_p1["SS_name_p1"]),
                                index=tmp_p1.index)
        tmp["Num_good_links_m1"].fillna(0, inplace=True)
        tmp["Num_good_links_p1"].fillna(0, inplace=True)
        tmp["Num_good_links"] = tmp["Num_good_links_m1"] + tmp["Num_good_links_p1"]
        
        # Add an assignment confidence column
        tmp["Conf_m1"] = "N"
        tmp.loc[tmp["Num_good_links_m1"]>0,"Conf_m1"] = "W"
        tmp.loc[tmp["Num_good_links_m1"]>1,"Conf_m1"] = "S"
        tmp.loc[tmp["Max_mismatch_m1"]>threshold,"Conf_m1"] = "X"
        
        tmp["Conf_p1"] = "N"
        tmp.loc[tmp["Num_good_links_p1"]>0,"Conf_p1"] = "W"
        tmp.loc[tmp["Num_good_links_p1"]>1,"Conf_p1"] = "S"
        tmp.loc[tmp["Max_mismatch_p1"]>threshold,"Conf_p1"] = "X"
        
        tmp["Conf2"] = tmp["Conf_m1"]+tmp["Conf_p1"]
        # Reorder letters eg WS -> SW
        tmp["Conf2"] = tmp["Conf2"].replace(["WS","NS","XS","NW","XW","XN"], 
                                           ["SW","SN","SX","WN","WX","NX"])

        tmp["Confidence"] = tmp["Conf2"]
        tmp["Confidence"] = tmp["Confidence"].replace(["SS","SW","SN","WW"],"High")
        tmp["Confidence"] = tmp["Confidence"].replace(["SX","WN"],"Medium")
        tmp["Confidence"] = tmp["Confidence"].replace("WX","Low")
        tmp["Confidence"] = tmp["Confidence"].replace(["NX","XX"],"Unreliable")
        tmp["Confidence"] = tmp["Confidence"].replace("NN","Undefined")
        
        return(tmp[["SS_name","Res_name","Max_mismatch_m1","Max_mismatch_p1","Max_mismatch",
                    "Num_good_links_m1","Num_good_links_p1","Num_good_links","Confidence"]])
            
    def add_consistency_info(self, assign_df=None, threshold=0.2):
        """ Find maximum mismatch and number of 'significant' mismatches for 
        each residue
        
        threshold: Minimum carbon shift difference for sequential residues to
            count as mismatched
        """
        
        # If the user hasn't specified an assign_df, use one already calculated 
        # for this SNAPS_assigner instance
        if assign_df is None:
            set_assign_df = True
            assign_df = self.assign_df
        else:
            set_assign_df = False
        
        # First check if there are any sequential atoms
        carbons = pd.Series(["C","CA","CB"])
        carbons_m1 = carbons + "_m1"
        seq_atoms = carbons[carbons.isin(assign_df.columns) & 
                            carbons_m1.isin(assign_df.columns)]
    
        if seq_atoms.size==0:
            # You can't do a comparison
            assign_df["Max_mismatch_m1"] = np.NaN
            assign_df["Max_mismatch_p1"] = np.NaN
            assign_df["Num_good_links_m1"] = 0
            assign_df["Num_good_links_p1"] = 0
            assign_df["Confidence"] = "Undefined"
        else:
            tmp = self.check_matching_consistency(assign_df, threshold)
            assign_df = assign_df.merge(tmp, how="left")
            
            
        if set_assign_df:
            self.assign_df = assign_df
        return(assign_df)
    
    def find_alt_assignments(self, N=1, by_ss=True, verbose=False):
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
        alt_assignments by rank
        """

        log_prob_matrix = self.log_prob_matrix
        best_matching = self.assign_df.loc[:,["SS_name","Res_name"]]
        best_matching.index = best_matching["SS_name"]
        best_matching.index.name = None
        alt_matching = None
        
        # Calculate sum probability for the best matching
        best_sum_prob = self.calc_overall_matching_prob(best_matching)
        
        # Calculate the value used to penalise the best match for each residue
        penalty = 2*log_prob_matrix.min().min()     
        logging.debug("Penalty value: %f", penalty)
        
        # Initialise DataFrame for storing alt_assignments
        alt_matching_all = best_matching.copy()
        alt_matching_all["Rank"] = 1
        alt_matching_all["Rel_prob"] = 0
                      
        for i in best_matching.index:   # Consider each spin system in turn
            ss = best_matching.loc[i, "SS_name"]
            res = best_matching.loc[i, "Res_name"]
            logging.debug("Finding alt assignments for original match %s - %s", ss, res)
            if verbose: print(ss, res)
            
            excluded = best_matching.loc[[i], :]
            
            for j in range(N):
                alt_matching = self.find_best_assignments(exc=excluded)
                                
                alt_matching["Rank"] = j+2
                alt_sum_prob = self.calc_overall_matching_prob(alt_matching)
                alt_matching["Rel_prob"] = alt_sum_prob - best_sum_prob
                               
                # Add the alt match for this ss or res to the results dataframe 
                # and also the excluded dataframe.
                if by_ss:
                    alt_matching_all = alt_matching_all.append(
                            alt_matching.loc[alt_matching["SS_name"]==ss, :], 
                            ignore_index=True)
                    res = alt_matching.loc[alt_matching["SS_name"]==ss, 
                                           "Res_name"].tolist()[0]
                    # The .tolist()[0] is to convert a single-item series into a string.
                else:
                    alt_matching_all = alt_matching_all.append(
                            alt_matching.loc[alt_matching["Res_name"]==res, :], 
                            ignore_index=True)
                    ss = alt_matching.loc[alt_matching["Res_name"]==res, 
                                          "SS_name"].tolist()[0]
                excluded = excluded.append(pd.DataFrame({"SS_name":[ss],"Res_name":[res]}), 
                                           ignore_index=True)
                   
        self.alt_assign_df = self.make_assign_df(alt_matching_all)
        if by_ss:
            self.alt_assign_df = self.alt_assign_df.sort_values(
                                                by=["SS_name", "Rank"])
        else:
            self.alt_assign_df = self.alt_assign_df.sort_values(
                                                by=["Res_name", "Rank"])
            
        return(self.alt_assign_df)
    
    def find_kbest_assignments(self, k, init_inc=None, init_exc=None, verbose=False):
        """ Find the k best overall assignments using the Murty algorithm.
        
        k: the number of assignments to find
        verbose: print information about progress
        
        This algorithm works by defining nodes. A node is a particular set of 
        constraints on the assignment, consisting of pairs that muct be included,
        and pairs that must be excluded. These constraints define a set of 
        potential assignments, one of which will be optimal. The clever part is 
        that this set of assignments can be further subdivided into the optimal 
        assignment plus a set of child nodes, each with additional constraints.
        The optimal solution can be found for each child node, allowing them to 
        be ranked. The highest ranking child can then be split up further, 
        allowing computational effort to be focused on the sets containing the 
        best solutions.
        For more details, see: 
        Murty, K. (1968). An Algorithm for Ranking all the Assignments in 
        Order of Increasing Cost. Operations Research, 16(3), 682-687
        """
        Node = namedtuple("Node", ["sum_log_prob","matching","inc","exc"])
        
        # Initial best matching (subject to initial constraints)
        best_matching = self.find_best_assignments(inc=init_inc, exc=init_exc)
        best_matching.index = best_matching["SS_name"]
        best_matching.index.name = None
        
        # Define lists to keep track of nodes
        ranked_nodes = SortedListWithKey(key=lambda n: n.sum_log_prob)
        unranked_nodes = SortedListWithKey(
                            [Node(self.calc_overall_matching_prob(best_matching),
                            best_matching, inc=init_inc, exc=init_exc)],
                            key=lambda n: n.sum_log_prob)
        
        while len(ranked_nodes)<k:
            # Set highest scoring unranked node as current_node
            current_node = unranked_nodes.pop()
            
            if verbose:
                s = str(len(ranked_nodes))+"\t"
                if current_node.inc is not None:
                    s = s + "inc:" +str(len(current_node.inc))+ "\t"
                if current_node.exc is not None:
                    s = s + "exc:"+ str(len(current_node.exc))
                print(s)
            
            # If the current node has forced included pairings, get a list of 
            # all parts of the matching that can vary.
            if current_node.inc is not None:
                matching_reduced = current_node.matching[
                                            ~current_node.matching["SS_name"].
                                            isin(current_node.inc["SS_name"])]
            else:
                matching_reduced = current_node.matching
            
            # Create child nodes and add them to the unranked_nodes list
            # -1 in range(matching_reduced.shape[0]-1) is because forcing 
            # inclusion of all but one residue actually results in getting back 
            # the original matching, and also results in an infinite loop of 
            # identical child nodes...
            for i in range(matching_reduced.shape[0]-1):
                #print("\t", i)
                #Construct dataframes of excluded and included pairs
                exc_i = matching_reduced.iloc[[i],:]
                if current_node.exc is not None:
                    exc_i = exc_i.append(current_node.exc, ignore_index=True)
                inc_i = matching_reduced.iloc[0:i,:]
                if current_node.inc is not None:
                    inc_i = inc_i.append(current_node.inc, ignore_index=True)
                inc_i = inc_i.drop_duplicates()
                # If there are no included pairs, it's better for inc_i to be None than an empty df
                if inc_i.shape[0]==0: 
                    inc_i = None
                
                matching_i = self.find_best_assignments(inc=inc_i, exc=exc_i, 
                                                     return_none_if_all_dummy=True,
                                                     verbose=False)
                if matching_i is None:
                    # If the non-constrained residues or spin systems are all 
                    # dummies, find_best_assignments will return None, and this 
                    # node can be discarded
                    pass
                else:
                    # Create a new child node and add to unranked_nodes
                    matching_i.index = matching_i["SS_name"]
                    matching_i.index.name = None
                    node_i = Node(self.calc_overall_matching_prob(matching_i), 
                                  matching_i, inc_i, exc_i)
                    unranked_nodes.add(node_i)
            ranked_nodes.add(current_node)
            
        return(ranked_nodes, unranked_nodes)                                                                                                                             
        
    def find_consistent_assignments(self, search_depth=10):
        """Find a consistent set of assignments using kbest search
        
        Finds the best assignment, then holds constant any consistent residues. Then finds 
        the k-best alternative assignments, and checks if a more consistent assignment 
        is found. If so, the extra consistent residues are held constant, and the 
        process is repeated.
        
        """
        
        matching0 = self.find_best_assignments()
        i=1
        
        while True:
            consistency0 = self.check_matching_consistency(matching0)
            consistency_sum0 = sum((consistency0["Max_mismatch"]<0.1) * 2**consistency0["Num_good_links"])
            print(i, consistency_sum0)
            i += 1
            
            mask = (consistency0["Max_mismatch"]<0.2) & (consistency0["Num_good_links"]>=4)
            inc0 = consistency0.loc[mask,["SS_name","Res_name"]]
        
            ranked, unranked = self.find_kbest_assignments(search_depth, init_inc=inc0, verbose=True)

            tmp = []
            for x in ranked:
                tmp2 = self.check_matching_consistency(x.matching)
                tmp += [sum((tmp2["Max_mismatch"]<0.2) * 2**tmp2["Num_good_links"])]
            
            if max(tmp) > consistency_sum0:
                # Prepare to loop back
                matching0 = ranked[tmp.index(max(tmp))].matching
                pass
            else:
                break
        
        sum(matching0["SS_name"] == matching0["SS_name"])
        
        return(matching0)
#        init_consistency = self.check_matching_consistency(new_matching)
#        old_consistency_sum = sum((init_consistency["Max_mismatch"]<0.1) * 2**init_consistency["Num_good_links"])
#        
#        mask = (init_consistency["Max_mismatch"]<0.1) & (init_consistency["Num_good_links"]>=4)
#        init_inc = init_consistency.loc[mask,["SS_name","Res_name"]]
#        
#        ranked, unranked = self.find_kbest_assignments(search_depth, init_inc=init_inc, verbose=True)
#        tmp = []
#        for x in ranked:
#            tmp2 = self.check_matching_consistency(x.matching)
#            tmp += [sum((tmp2["Max_mismatch"]<0.1) * 2**tmp2["Num_good_links"])]
#        
#        if max(tmp) > old_consistency_sum:
#            pass
#        else:
#            break
#        new_matching = ranked[tmp.index(max(tmp))].matching
    
    def find_consistent_assignments2(self, threshold=0.2):
        """Try to find a consistent set of assignments by optimising both match 
        to predictions and mismatches between adjacent residues"""
        matching0 = self.find_best_assignments()
        assign_df0 = self.make_assign_df(matching0)
        assign_df0 = self.add_consistency_info(assign_df0, threshold)
        N_HM_conf0 = assign_df0["Confidence"].isin(["High","Medium"]).sum()
        print(N_HM_conf0)
        
        while True:
            # Find all residues with an adjacent confident residue
            HM_conf_res = assign_df0.loc[assign_df0["Confidence"].isin(["High","Medium"]),
                                         "Res_name"]
            
            #Make a dataframe for looking up which spin system is assigned to each residue
            tmp = matching0
            tmp.index = tmp["Res_name"]
            
            # Limit their assignment options to consistent spin systems
            a = deepcopy(self)
            
            for res in HM_conf_res:
                # Get the neighbouring residues
                res_m1 = self.neighbour_df.loc[res,"Res_name_m1"]
                res_p1 = self.neighbour_df.loc[res,"Res_name_p1"]
                
                ss = tmp.SS_name[res]
                
                # (Need to be careful with NaN values at this point)
                if not pd.isna(res_m1):
                    #Get list of inconsistent i-1 spin systems
                    disallowed_ss = self.mismatch_matrix.index[
                            self.mismatch_matrix.loc[:,ss]>threshold]
                    # Set the inconsistent spins systems to have a high log_prob
                    a.log_prob_matrix.loc[disallowed_ss, res_m1] = a.log_prob_matrix.min().min()
                if not pd.isna(res_p1):
                    #Get list of inconsistent i+1 spin systems
                    disallowed_ss = self.mismatch_matrix.columns[
                            self.mismatch_matrix.loc[ss,:]>threshold]
                    #Set the inconsistent spins systems to have a high log_prob
                    a.log_prob_matrix.loc[disallowed_ss, res_p1] = a.log_prob_matrix.min().min()


            # Rerun assignment and see how many are consistent now
            #return(a.log_prob_matrix)
            matching1 = a.find_best_assignments()
            assign_df1 = a.make_assign_df(matching1)
            assign_df1 = a.add_consistency_info(assign_df1, threshold)
            N_HM_conf1 = assign_df1["Confidence"].isin(["High","Medium"]).sum()
            print(N_HM_conf1)
            
            if N_HM_conf1 > N_HM_conf0:
                assign_df0 = assign_df1
                N_HM_conf0 = N_HM_conf1
            else:
                break

        return(assign_df0)
    
    def find_seq_assignment(self):
        """Find the ordering that maximises the number of good sequential links"""
        
        row_ind, col_ind = linear_sum_assignment(self.mismatch_matrix-1*self.consistent_links_matrix)
        matching = pd.DataFrame({"i_m1":self.consistent_links_matrix.index[row_ind],
                                 "i":self.consistent_links_matrix.columns[col_ind]})

        
        return(matching)
    
    def output_shiftlist(self, filepath, format="sparky", 
                         confidence_list=["High","Medium","Low","Unreliable","Undefined"]):
        """Export a chemical shift list for use with other programs
        """
        atoms = {"H","N","HA","C","CA","CB"}.intersection(self.assign_df.columns)
        
        df_wide = self.assign_df.loc[self.assign_df["Confidence"].isin(confidence_list),["Res_N","Res_type"]+list(atoms)]
        # Check in case no shifts were selected
        if df_wide.empty:   
            print("No chemical shifts selected for export.")
            # Create an empty file
            df_wide.to_csv(filepath, sep="\t", float_format="%.3f",
                           index=True, header=False)
            return(None)
        
        df = df_wide.melt(id_vars=["Res_N","Res_type"], value_vars=atoms, 
                          var_name="Atom_type", value_name="Shift")
        df = df.sort_values(["Res_N","Atom_type"])
        
        if format=="sparky":
            df["Group"] = df["Res_type"] + df["Res_N"].astype(str)
            
            df = df.rename(columns={"Atom_type":"Atom"})
            df.loc[df["Atom"]=="H","Atom"] = "HN"
            
            nuc_dict = {"HN":"1H", "N":"15N", "HA":"1H",
                        "C":"13C", "CA":"13C", "CB":"13C"}
            df["Nuc"] = [nuc_dict[a] for a in df["Atom"]]
            
            output_df = df[["Group","Atom","Nuc","Shift"]]
            output_df["Sdev"] = 0.0
            output_df["Assignments"] = int(1)
            
            output_df = output_df.dropna()
            
            if filepath is not None:
                output_df.to_csv(filepath, sep="\t", float_format="%.3f",
                                 index=False)
        elif format=="xeasy":
            df.loc[df["Atom_type"]=="H","Atom_type"] = "HN"
            
            output_df = df[["Shift","Atom_type","Res_N"]]
            output_df.insert(1, "Sdev", 0)
            
            output_df["Shift"] = output_df["Shift"].fillna(999.0)
            output_df = output_df.dropna()
            output_df = output_df.reset_index(drop=True)
        
            if filepath is not None:
                output_df.to_csv(filepath, sep="\t", float_format="%.3f",
                                 index=True, header=False)
        elif format=="nmrpipe":
            df.loc[df["Atom_type"]=="H","Atom_type"] = "HN"
            
            output_df = df[["Res_N", "Res_type", "Atom_type", "Shift"]]
            
            # Prepare the sequence info
            #TODO: This doesn't work properly if some residues are missing 
            # (eg. if only high confidence assignments are being exported). 
            tmp = pd.DataFrame({"Res_N":np.arange(df["Res_N"].min(), 
                                                  df["Res_N"].max()+1)})
            tmp = tmp.merge(df_wide[["Res_N","Res_type"]], how="left", on="Res_N")
            tmp["Res_type"] = tmp["Res_type"].fillna("X")
            seq = tmp["Res_type"].sum()
            
            if filepath is not None:
                f = open(filepath, 'w+')
                f.write("REMARK Chemical shifts (automatically assigned using SNAPS)\n\n")
                f.write("DATA FIRST_RESID %d\n\n" % tmp["Res_N"].min())
                f.write("DATA SEQUENCE %s\n\n" % seq)
                f.write("VARS   RESID RESNAME ATOMNAME SHIFT\n")
                f.write("FORMAT %4d   %1s     %4s      %8.3f\n\n")
                
                for i in output_df.index:
                    f.write("%4d %1s %4s %8.3f\n" % (output_df.loc[i, "Res_N"],
                                                     output_df.loc[i, "Res_type"],
                                                     output_df.loc[i, "Atom_type"],
                                                     output_df.loc[i, "Shift"]))
                f.close()
                
        else:
            print("format string '%s' no recognised." % format)
            return(None)
            
        return(output_df)
    
    def output_peaklists(self, filepath, format="sparky", 
                         spectra=["hsqc","hnco","hncaco","hncacb", "hncocacb"]):
        """ Output assigned peaklists
        """
        return(0)
    
    
    def plot_strips(self, outfile=None, format="html", return_json=True, plot_width=1000):
        """Make a strip plot of the assignment.
        
        Uses bokeh module for plotting. Returns the bokeh plot object.
        
        outfile: if defined, the plot will be saved to this location
        format: type of output. Either "html" or "png"
        """
        df = self.assign_df
        plotlist = []
        
        # First check if there are any sequential atoms
        carbons = pd.Series(["C","CA","CB"])
        carbons_m1 = carbons + "_m1"
        seq_atoms = carbons[carbons.isin(df.columns) & 
                            carbons_m1.isin(df.columns)]
        
        if seq_atoms.size==0:
            # You can't draw a strip plot
            print("No sequential links in data - strip plot not drawn.")
            return(None)
        else:
            #Add extra rows for missing residues (eg Prolines)
            tmp = list(range(int(df["Res_N"].min()), int(df["Res_N"].max()+1)))
            tmp2 = pd.DataFrame({"Res_N":tmp, 
                                 "Dummy_res":[False]*len(tmp)})
            df = pd.merge(tmp2,df, on=None, how="outer")
            
            mask = df["Res_name"].isna()
            df.loc[mask,"Res_name"] = df.loc[mask,"Res_N"].astype(int).astype(str) + "_"
            df.loc[mask,"Dummy_SS"] = True
            df.loc[mask,"Dummy_res"] = False
            #df.loc[mask,"SS_name"] = df.loc[mask,"Res_N"].astype(str) + "_"
            df.loc[mask,"Res_name"] = [s.rjust(5) for s in df.loc[mask,"Res_name"]]
            
            df["Dummy_res"] = df["Dummy_res"].astype(bool)
            
            tmp_plt = figure(x_range=df["Res_name"].tolist())
            for atom in seq_atoms:
                # Setup plot
                plt = figure(y_axis_label=atom+" (ppm)",
                             x_range=tmp_plt.x_range,
                             tools="xpan, xwheel_zoom,save,reset",
                             height=200, width=plot_width)
                #plt.toolbar.active_scroll = plt.select_one(WheelZoomTool) 
                
                ## Plot the vertical lines
                vlines = df.loc[~df["Dummy_res"], ["Res_name",atom,atom+"_m1"]]
                # Introduce NaN's to break the line into discontinuous segments
                # "_z" in name is to ensure it sorts after the atom+"_m1" shifts
                # eg. CA, CA_m1, CA_z
                vlines[atom+"_z"] = np.NaN  
                # Convert from wide to long
                vlines = vlines.melt(id_vars=["Res_name"], 
                                     value_vars=[atom, atom+"_m1", atom+"_z"], 
                                     var_name="Atom_type", 
                                     value_name="Shift")
                vlines = vlines.sort_values(["Res_name","Atom_type"], 
                                            ascending=[True,False])
                plt.line(vlines["Res_name"], vlines["Shift"], 
                         line_color="black", line_dash="dashed")
                
                ## Plot the horizontal lines
                hlines = df.loc[~df["Dummy_res"], ["Res_N","Res_name",atom,atom+"_m1"]]
                # Introduce NaN's to break the line into discontinuous segments
                # "_a" in name ensures it sorts between the atom and atom+"_m1" shifts
                # eg. CA, CA_a, CA_m1
                hlines[atom+"_a"] = np.NaN  
                # Add extra lines for missing residues
#                tmp = list(range(hlines["Res_N"].min(), hlines["Res_N"].max()))
#                tmp2 = pd.DataFrame({"Res_N":tmp})
#                hlines = pd.merge(tmp2,hlines, on=None)
#                mask = hlines["Res_name"].isna()
#                hlines.loc[mask,"Res_name"] = hlines.loc[mask,"Res_N"].astype(str)
#                hlines["Res_name"] = [s.rjust(5) for s in hlines["Res_name"]]
#                # Convert from wide to long
                hlines = hlines.melt(id_vars=["Res_name"], 
                                     value_vars=[atom, atom+"_m1", atom+"_a"], 
                                     var_name="Atom_type", 
                                     value_name="Shift")
                hlines = hlines.sort_values(["Res_name","Atom_type"], 
                                            ascending=[True,False])
                plt.line(hlines["Res_name"], hlines["Shift"], 
                         line_color="black", line_dash="solid")
                
                # Draw circles at the observed chemical shifts
                plt.circle(df["Res_name"], df[atom], fill_color="blue", size=5)
                plt.circle(df["Res_name"], df[atom+"_m1"], fill_color="red", size=5)
                
                # Reverse the y axis and set range:
                plt.y_range = Range1d(
                        df[[atom,atom+"_m1"]].max().max()+5, 
                        df[[atom,atom+"_m1"]].min().min()-5)
                
                # Change axis label orientation
                plt.xaxis.major_label_orientation = 3.14159/2
                plt.xaxis.visible = False
                
                plotlist = plotlist + [plt]
            
            plotlist[-1].xaxis.visible = True
            
            #### Make mismatch plot ####
            plt = figure(title="Mismatch plot",
                                y_axis_label="Mismatch (ppm)",
                                x_range=tmp_plt.x_range,
                                tools="xpan, xwheel_zoom,save,reset",
                                height=200, width=plot_width)
            #plt.toolbar.active_scroll = plt.select_one(WheelZoomTool) 
            
            plt.vbar(x=df["Res_name"], 
                     top=df[["Max_mismatch_m1","Max_mismatch_p1"]].max(axis=1), 
                     width=1)
            
            # Draw a line showing the threshold for mismatches
            threshold_line = Span(location=self.pars["seq_link_threshold"], 
                                  dimension="width", line_dash="dashed", 
                                  line_color="red")
            plt.add_layout(threshold_line)
            
            #plt.vbar(x=df["Res_name"], top=-df["Max_mismatch_m1"], width=1)
            # Change axis label orientation
            plt.xaxis.major_label_orientation = 3.14159/2
            plotlist = [plt] + plotlist
            
            # Make confidence plot
            plt = figure(title="Confidence plot (mouseover to see observed spin system name)",
                                x_range=tmp_plt.x_range,
                                tools="xpan, xwheel_zoom,hover,save,reset",
                                tooltips=[("Pred", "@Res_name"),("Obs","@SS_name")],
                                height=100, width=plot_width)
            #plt.toolbar.active_scroll = plt.select_one(WheelZoomTool)
            
            # Create a colour map based on confidence
            colourmap = {"High":"green",
                         "Medium":"yellowgreen",
                         "Low":"orange",
                         "Unreliable":"red",
                         "Undefined":"grey"}
            
            # Plot the peaks
            for k in colourmap.keys():
                tmp = ColumnDataSource(df[df["Confidence"]==k])
                plt.vbar(x="Res_name", top=1, width=1, 
                         color=colourmap[k], legend=k, source=tmp)
            
            plt.legend.orientation = "horizontal"
            plt.legend.location = "top_center"
            # Change axis label orientation
            #plt.xaxis.major_label_orientation = 3.14159/2
            plt.axis.visible = False
            
            plotlist = [plt] + plotlist
            
            # Put them all together, and export
            p = gridplot(plotlist, ncols=1)

            if outfile is not None:
                Path(outfile).resolve().parents[1].mkdir(parents=True, exist_ok=True)
                if format=="html":
                    output_file(outfile)
                    save(p)
                elif format=="png":
                    export_png(p, outfile)
            
            if return_json:
                return(json_item(p, "strip_plot"))
            else:
                return(p)
                
    def plot_hsqc(self, outfile=None, format="html", return_json=True):
        """Plot an HSQC coloured by prediction confidence"""
        
        assign_df = self.assign_df.copy()
        
        plt = figure(title="HSQC",
                        x_axis_label="1H (ppm)",
                        y_axis_label="15N (ppm)",
                        height=750, width=750)
        
        # Create a colour map based on confidence
        colourmap = {"High":"green",
                         "Medium":"yellowgreen",
                         "Low":"orange",
                         "Unreliable":"red",
                         "Undefined":"grey"}

        
        # Plot the peaks
        for k in colourmap.keys():
            tmp = assign_df[assign_df["Confidence"]==k]
            plt.circle(tmp["H"], tmp["N"], color=colourmap[k], legend=k)
        
        # Label the points
        df = ColumnDataSource(assign_df)
        labels = LabelSet(x="H", y="N", text="Res_name", source=df, text_font_size="10pt")
        plt.add_layout(labels)
        
        # Reverse the axes and set range
        plt.x_range = Range1d(assign_df["H"].max()+1, assign_df["H"].min()-1)
        plt.y_range = Range1d(assign_df["N"].max()+5, assign_df["N"].min()-5)
        
        # Change legend layout
        plt.legend.orientation = "horizontal"
        
        if outfile is not None:
            Path(outfile).resolve().parents[1].mkdir(parents=True, exist_ok=True)
            if format=="html":
                output_file(outfile)
                save(plt)
            elif format=="png":
                export_png(plt, outfile)
        
        if return_json:
            return(json_item(plt))
        else:
            return(plt)
    
    def plot_strips_plotnine(self, atom_list=["C","C_m1","CA","CA_m1","CB","CB_m1"]):
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
        plot_df.loc[plot_df["Atom_type"].isin(["C_m1","CA_m1","CB_m1"]),"i"] = "-1"
        plot_df["Atom_type"] = plot_df["Atom_type"].replace({"C_m1":"C", 
                                                   "CA_m1":"CA", "CB_m1":"CB"}) 
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
        if not all(pd.Series(['Max_mismatch_m1', 'Max_mismatch_p1']).
                   isin(assign_df.columns)):
            return(None)
        else:
            # Pad Res_name column with spaces so that sorting works correctly
            assign_df["Res_name"] = assign_df["Res_name"].str.pad(6)
            assign_df["x_name"] = (assign_df["Res_name"] + "_(" + 
                                     assign_df["SS_name"] + ")")
            
            # Make the plot
            plt = ggplot(aes(x="x_name"), data=assign_df) 
            plt = plt + geom_col(aes(y="abs(Max_mismatch_m1)"))
            plt = plt + xlab("Residue name")
            plt = plt + ylab("Mismatch to previous residue (ppm)")
            plt = plt + theme_bw() + theme(axis_text_x = element_text(angle=90))
                   
            return(plt)
        