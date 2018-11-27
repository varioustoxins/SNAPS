#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 16:12:52 2018

@author: aph516
"""

import numpy as np
import pandas as pd
from scipy.stats import norm, multivariate_normal
from plotnine import *
from NAPS_assigner import NAPS_assigner
from pathlib import Path
#path = "/Users/aph516/GitHub/NAPS/"
path = Path("C:/Users/Alex/GitHub/NAPS/")

#%%
#### Prepare to import all proteins in the testset

testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/"data/testset/simplified_BMRB"/file for file in testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/"data/testset/shiftx2_results"/file for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df.index = testset_df["ID"]

a = NAPS_assigner()

# Work out which residue types don't match between obs and preds
seq_df = None
for i in testset_df["ID"]:
    obs = a.import_obs_shifts(testset_df.loc[i, "obs_file"])
    # Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    
    # Make columns for the i-1 predicted shifts of C, CA and CB
    obs.index = obs["Res_N"]
    obs_m1 = obs[["Res_type","Res_N"]].copy()
    obs_m1.index = obs_m1.index+1
    obs_m1.columns = obs_m1.columns + "m1"
    obs = pd.merge(obs, obs_m1, how="left", left_index=True, right_index=True)
    obs.index = obs["SS_name"]
    
    
    preds = a.read_shiftx2(testset_df.loc[i, "preds_file"])
    #Change Bs to Cs
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_typem1"]=="B", "Res_typem1"] = "C"
    
    df = pd.merge(obs, preds, on=["Res_N"], how="outer", suffixes=["_obs","_pred"])
    df["ID"] = i
    if type(seq_df)!=type(df):
        seq_df = df.copy()
    else:
        seq_df = pd.concat([seq_df, df], ignore_index=True)

seq_df = seq_df[["ID","Res_N","Res_type_obs","Res_type_pred","Res_typem1_obs","Res_typem1_pred"]]
seq_df = seq_df.fillna("NA")
tmp1 = seq_df.loc[seq_df["Res_type_obs"]!=seq_df["Res_type_pred"],:]
tmp1 = tmp1.loc[~tmp1["Res_type_obs"].isin(["NA"]) & ~tmp1["Res_type_pred"].isin(["NA"]),:]
tmp2 = seq_df.loc[seq_df["Res_typem1_obs"]!=seq_df["Res_typem1_pred"],:]
tmp2 = tmp2.loc[~tmp2["Res_typem1_obs"].isin(["NA"]) & ~tmp2["Res_typem1_pred"].isin(["NA"]),:]

bad_res = pd.concat([tmp1, tmp2])
bad_res = bad_res.drop_duplicates()

#%%

#### Import the actual data
obs_all = None
for i in testset_df["ID"]:
#for i in ["A063"]:
    obs = a.import_obs_shifts(testset_df.loc[i, "obs_file"])
    preds = a.read_shiftx2(testset_df.loc[i, "preds_file"])
    
    #Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_typem1"]=="B", "Res_typem1"] = "C"
    
    obs["ID"] = i
    preds["ID"] = i

    if type(obs_all)!=type(obs):
        obs_all = obs.copy()
        preds_all = preds.copy()
    else:
        obs_all = pd.concat([obs_all, obs], ignore_index=True)
        preds_all = pd.concat([preds_all, preds], ignore_index=True)

df = pd.merge(obs_all, preds_all, on=["ID","Res_N","Res_type"], how="outer", suffixes=["_obs","_pred"])

# Get rid of any residues with mismatching types, or where the type is NaN
df = df.loc[~(df["ID"]+df["Res_N"].astype(str)).isin(bad_res["ID"]+bad_res["Res_N"].astype(str)),:]
df = df.dropna(axis=0, subset=["Res_type","Res_typem1"])
# Make a normalised version of observed shift
# Do this by subtracting the mean shift, separately for each residue type
atoms_obs = [s+"_obs" for s in ["C","CA","CB","H","HA","N"]]
atoms_pred = [s+"_pred" for s in ["C","CA","CB","H","HA","N"]]
atoms_m1_obs = [s+"m1_obs" for s in ["C","CA","CB","H","HA","N"]]
atoms_m1_pred = [s+"m1_pred" for s in ["C","CA","CB","H","HA","N"]]
for r in df["Res_type"].unique():
    mu = df.loc[df["Res_type"]==r, atoms_obs].mean()
    sig = df.loc[df["Res_type"]==r, atoms_obs].std()
    # i shifts
    df.loc[df["Res_type"]==r, atoms_obs] = (df.loc[df["Res_type"]==r, atoms_obs] - mu)/sig
    mu.index = atoms_pred
    sig.index = atoms_pred
    df.loc[df["Res_type"]==r, atoms_pred] = (df.loc[df["Res_type"]==r, atoms_pred] - mu)/sig
    # i-1 shifts
    mu.index = atoms_m1_obs
    sig.index = atoms_m1_obs
    df.loc[df["Res_typem1"]==r, atoms_m1_obs[0:3]] = (df.loc[df["Res_typem1"]==r, atoms_m1_obs[0:3]] - mu[atoms_m1_obs[0:3]])/sig[atoms_m1_obs[0:3]]
    mu.index = atoms_m1_pred
    sig.index = atoms_m1_pred
    df.loc[df["Res_typem1"]==r,atoms_m1_pred[0:3]] = (df.loc[df["Res_typem1"]==r, atoms_m1_pred[0:3]] - mu[atoms_m1_pred[0:3]])/sig[atoms_m1_pred[0:3]]
        

# Calculate the prediction errors
for atom in a.pars["atom_set"]:
    df["d_"+atom] = df[atom+"_pred"] - df[atom+"_obs"]

#%%
#### Try out Yeo-Johnson transformation
# This should transform a distribution to be more normal

from YeoJohnson import YeoJohnson

yj = YeoJohnson()
test_data = df["d_CA"]



# To work out the best value of lambda, we'll need to find the maximum likelihood value of lambda.
results = pd.DataFrame(columns=["lam", "log_prob"], index=np.linspace(0,2,50))
for l in np.linspace(0,2,50):
    x = yj.fit(df["C_obs"].dropna(), l)
    mu, sig = norm.fit(x)
    p = sum(np.log(norm.pdf(x, mu, sig)))
    results.loc[l,:] = [l, p]
ggplot(data=results) + geom_line(aes(x="lam", y="log_prob", group=1))
    
trans_data = pd.DataFrame( {"test":df["C_obs"], "trans" : yj.fit(df["C_obs"], 1.1) })
trans_data = pd.DataFrame( {"test":df["d_C"], "trans" : yj.fit(df["d_C"], 1.1) })

ggplot(data=trans_data) + geom_density(aes(x="test")) + geom_density(aes(x="trans"), colour="red")


#%%
#### Try out multivariate normal distribution
df2 = df[["ID","SS_name","Res_N","Res_type","Res_typem1"]+[s+"_obs" for s in a.pars["atom_set"]]+["d_"+s for s in a.pars["atom_set"]]]

b = abs(df2.corr())
b[b<0.1]=-1

df2_mean = df2.mean()
df2_mean = df2_mean.drop("Res_N")
df2_cov = df2.cov()
df2_cov = df2_cov.drop(index="Res_N", columns="Res_N")
#df2_cov[abs(df2_cov)<0.1] = 0 # Can't set small covariances to zero because it makes matrix non-invertible, which causes error in multivariate_normal()


mvn_full = multivariate_normal(df2_mean, df2_cov)
mvn_obs = multivariate_normal(df2_mean.iloc[0:9], df2_cov.iloc[0:9, 0:9])   # The expected distribution of the observations
mvn_d = multivariate_normal(df2_mean.iloc[9:18], df2_cov.iloc[9:18, 9:18])  # This is the marginal distribution of the errors over the observations

# Test random variables drawn from distribution
tmp = pd.DataFrame(mvn_full.rvs(size=10000), columns = df2_cov.columns)
ggplot(data=tmp) + geom_point(aes(x="d_CBm1", y="CBm1_obs"))    # Distribution is more compressed than real one. I guess this is because the real distribution has longer tails than the normal
c = tmp.cov()
d = df2_cov - c     # The covariance of the results matches the original data fairly well at least

# How to add in conditional probability?
obs1 = [0]*18
mvn_full.pdf(obs1)

mvn_d.pdf(obs1[9:18])

from math import exp
e = exp(mvn_full.logpdf(obs1) - mvn_obs.logpdf(obs1[0:9]))


#%%
#### Exploratory analysis
ggplot(data=df2) + geom_point(aes(x="d_CBm1", y="CBm1_obs", colour="Res_typem1"))# + facet_wrap("Res_type")
#ggplot(data=df[df["Res_type"]!="C"])

tmp = df2.loc[df2["CBm1_obs"]>7.5,:]


leu = df.loc[df["Res_type"]=="L", ["d_"+atom for atom in a.pars["atom_set"]]]

#ggplot(data=leu) + geom_point(aes(x="CB_obs", y="CB_pred"))

leu_mean = leu.mean(axis=0)
leu_cov = leu.cov()
leu_corr = leu.corr()

delta = df.loc[:, ["ID","Res_name","Res_type", "Res_typem1"]+["d_"+atom for atom in a.pars["atom_set"]]]
d_mean = delta.mean(axis=0)
d_cov = delta.cov()
d_corr = delta.corr()

ggplot(data=delta) + geom_point(aes(x="d_CAm1", y="d_CBm1", colour="Res_type"))



tmp = delta.loc[abs(delta["d_CBm1"])>5,:]
tmp2 = df.loc[(df["ID"]+df["Res_name"]).isin(tmp["ID"]+tmp["Res_name"]),:]

