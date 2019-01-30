#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:51:28 2019

@author: aph516
"""

import numpy as np
import pandas as pd
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
from pathlib import Path
from distutils.util import strtobool
from scipy.stats import norm
from math import isnan, log10

path = Path("/Users/aph516/GitHub/NAPS/")

a = NAPS_assigner()
    
# Import config file
config = pd.read_table(path/"config/config.txt", sep="\s+", comment="#", header=None,
                       index_col=0, names=["Value"])
a.pars["pred_offset"] = int(config.loc["pred_offset"].Value)
a.pars["prob_method"] = config.loc["prob_method"].Value
a.pars["alt_assignments"] = int(config.loc["alt_assignments"].Value)
a.pars["atom_set"] = {s.strip() for s in config.loc["atom_set"].Value.split(",")}
tmp = [s.strip() for s in config.loc["atom_sd"].Value.split(",")]
a.pars["atom_sd"] = dict([(x.split(":")[0], float(x.split(":")[1])) for x in tmp])
a.pars["plot_strips"] = bool(strtobool(config.loc["plot_strips"].Value))


# Import observed and predicted shifts
importer = NAPS_importer()


importer.import_testset_shifts(path/"data/testset/simplified_BMRB/4032.txt")

a.obs = importer.obs

a.import_pred_shifts(path/"data/testset/shiftx2_results/A001_1KF3A.cs", "shiftx2")

# Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix(sf=1, verbose=False)
assign_df, best_match_indexes = a.find_best_assignment()
a.check_assignment_consistency(threshold=0.1)

#%% Test stuff

obs = a.obs
preds = a.preds

obs_H = obs["H"].repeat(len(obs.index)).values.reshape([len(obs.index),-1])
preds_H = preds["H"].repeat(len(preds.index)).values.reshape([len(preds.index),-1]).transpose()
delta_H = preds_H - obs_H
prob_H = norm.logpdf(delta_H.to_numeric())

#%%
def calc_log_prob_matrix2(assigner, default_prob=0.01):
    # Use default atom_sd values if not defined
    atom_sd = assigner.pars["atom_sd"]
    
    atoms = assigner.pars["atom_set"].intersection(obs.columns)
    
    log_prob_matrix = pd.DataFrame(0, index=obs.index, columns=preds.index)
    
    for atom in atoms:
        obs_atom = obs[atom].repeat(len(obs.index)).values.reshape([len(obs.index),-1])
        preds_atom = preds[atom].repeat(len(preds.index)).values.reshape([len(preds.index),-1]).transpose()
        
        #return(obs[atom])
        delta_atom = preds_atom - obs_atom
        
        # Make a note of NA positions in delta, and set them to zero 
        # (this avoids warnings when using norm.cdf later)
        na_mask = np.isnan(delta_atom)
        delta_atom[na_mask] = 0
        
        prob_atom = pd.DataFrame(norm.logpdf(delta_atom, scale=atom_sd[atom]),
                                 index=obs.index, columns=preds.index)
        prob_atom[na_mask] = log10(default_prob)
        
        log_prob_matrix = log_prob_matrix + prob_atom
    
    return(log_prob_matrix)

a.pars["prob_method"]="cdf"
lpm = a.log_prob_matrix
lpm2 = calc_log_prob_matrix2(a)
lpm3 = a.calc_log_prob_matrix2()
