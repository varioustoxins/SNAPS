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
id = "A069"

importer1 = NAPS_importer()
importer1.import_testset_shifts(testset_df.loc[id, "obs_file"])

importer2 = NAPS_importer()
#importer2.import_testset_shifts(testset_df.loc[id, "obs_file"])
AA_class, AA_classm1 = "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML".split(";")
importer2.import_testset_shifts(testset_df.loc[id, "obs_file"],
                               SS_class=AA_class.split(","),
                               SS_classm1=AA_classm1.split(","))

#%%
a = NAPS_assigner()
b = NAPS_assigner()

a.obs = importer1.obs
b.obs = importer2.obs

# Import config file
a.read_config_file(path/"config/config.txt")
b.read_config_file(path/"config/config_hadamac.txt")

a.import_pred_shifts(testset_df.loc[id, "preds_file"], "shiftx2")
b.import_pred_shifts(testset_df.loc[id, "noshifty_file"], "shiftx2")

# Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix2(sf=1, verbose=False)
matching = a.find_best_assignments()
assign_df = a.make_assign_df(matching, set_assign_df=True)
assign_df = a.check_assignment_consistency(threshold=0.1)
#alt_assign_df = a.find_alt_assignments(N=2, verbose=True)

b.add_dummy_rows()
tmp = b.calc_log_prob_matrix2(sf=1, verbose=False)
matching2 = b.find_best_assignments()
assign_df2 = b.make_assign_df(matching2, set_assign_df=True)
assign_df2 = b.check_assignment_consistency(threshold=0.1)
#alt_assign_df2 = b.find_alt_assignments(N=2, by_ss=True, verbose=True)

obs = a.obs
preds = a.preds
log_prob_matrix = a.log_prob_matrix


obs2 = b.obs
preds2 = b.preds
log_prob_matrix2 = b.log_prob_matrix
#preds_corr = b.preds_corr
#pc_CA = preds_corr["CA"]
#del_CA = pc_CA - pd.DataFrame(preds["CA"].repeat(len(preds.index)).values.
#                                reshape([len(preds.index),-1]).transpose(),
#                                index=obs.index, columns=preds.index)

sum(assign_df["SS_name"] == assign_df["Res_name"])
sum(assign_df2["SS_name"] == assign_df2["Res_name"])

tmp1 = assign_df
tmp1.index = tmp1["Res_name"]
tmp2 = assign_df2
tmp2.index = tmp2["Res_name"]

sum(tmp1["SS_name"] == tmp2["SS_name"])

tmp = log_prob_matrix.melt()

# Test kbest assignments
kbest, unranked = a.find_kbest_assignments(10, verbose=True)

tmp0 = kbest[0]
tmp1 = kbest[1]
tmp9 = kbest[9]
#[n.sum_log_prob for n in kbest]
#[n.sum_log_prob for n in unranked]
unranked[-1].exc
#tmp = alt_assign_df[["SS_name","Res_name","Rank","Rel_prob", "Log_prob"]]
#plt = a.plot_strips()
#%% Test stuff

obs = a.obs
preds = a.preds

obs_H = obs["H"].repeat(len(obs.index)).values.reshape([len(obs.index),-1])
preds_H = preds["H"].repeat(len(preds.index)).values.reshape([len(preds.index),-1]).transpose()
delta_H = preds_H - obs_H
prob_H = norm.logpdf(delta_H.to_numeric())

dist_mat = a.calc_dist_matrix(use_atoms=["H","N"], rank=True)
rank_dist = dist_mat.rank(axis=1)

common_residues = list(rank_dist.index[rank_dist.index.isin(rank_dist.columns)])
tmp = rank_dist.lookup(common_residues, common_residues)
#tmp = pd.Series(np.diag(rank_dist), index=rank_dist.index)

#%% Write a file with HADAMAC info
tmp = obs.loc[:,["SS_name", "SS_classm1"]]
tmp["Type"] = "in"
tmp = tmp.dropna()
tmp.to_csv(path/"data/SS_class_info.txt", sep="\t", header=False, index=False)
#%%
x=np.indices((3,3))
x=np.swapaxes(x,0,2)
from scipy.stats import multivariate_normal
mvn = multivariate_normal([0,0],[1,1])
mvn.pdf(x)

x=np.array([[0,1],[2,3]])
y=np.array([[4,5],[6,7]])
z=np.array([[8,9],[10,np.NaN]])
xyz = np.array([x,y,z])
xyz = np.moveaxis(xyz, 0, -1)
mvn = multivariate_normal([0,0,0],[1,1,1])
mvn.pdf(xyz)




#x2 = pd.DataFrame(x)
#y2 = pd.DataFrame(y)
#z2 = pd.DataFrame(z)
#
#xyz2 = np.array([x2,y2,z2])
