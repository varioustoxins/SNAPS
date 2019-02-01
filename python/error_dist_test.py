#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 16:12:52 2018

@author: aph516
"""

import numpy as np
import pandas as pd
from scipy.stats import norm, multivariate_normal, linregress
from plotnine import *
from plotnine.ggplot import save_as_pdf_pages
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
from pathlib import Path
from math import sqrt
path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/Alex/GitHub/NAPS/")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")
#%%
#### Prepare to import all proteins in the testset

testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/"data/testset/simplified_BMRB"/file for file in testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/"data/testset/shiftx2_results"/file for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df.index = testset_df["ID"]

importer = NAPS_importer()
assigner = NAPS_assigner()

# Work out which residue types don't match between obs and preds
seq_df = None
for i in testset_df["ID"]:
    obs = importer.import_testset_shifts(testset_df.loc[i, "obs_file"])
    # Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    
    # Make columns for the i-1 predicted shifts of C, CA and CB
    obs.index = obs["Res_N"]
    obs_m1 = obs[["Res_type","Res_N"]].copy()
    obs_m1.index = obs_m1.index+1
    obs_m1.columns = obs_m1.columns + "m1"
    obs = pd.merge(obs, obs_m1, how="left", left_index=True, right_index=True)
    obs.index = obs["SS_name"]
    
    preds = assigner.import_pred_shifts(testset_df.loc[i, "preds_file"], filetype="shiftx2")
    #Change Bs to Cs
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_typem1"]=="B", "Res_typem1"] = "C"
       
    df = pd.merge(obs, preds, on="Res_N", how="outer", suffixes=["_obs","_pred"])
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

importer = NAPS_importer()
assigner = NAPS_assigner()

atom_set = {"H","N","HA","C","CA","CB","Cm1","CAm1","CBm1"}

#### Import the actual data
obs_all = None
for i in testset_df["ID"]:
#for i in ["A063"]:
    obs = importer.import_testset_shifts(testset_df.loc[i, "obs_file"], 
                                         remove_Pro=False)
    preds = assigner.import_pred_shifts(testset_df.loc[i, "preds_file"], 
                                        filetype="shiftx2")
    
    #Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type"]=="B", "Res_name"] = (
            preds.loc[preds["Res_type"]=="B", "Res_name"].str.replace("B","C"))
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_typem1"]=="B", "Res_typem1"] = "C"
    

    #Convert wide to long
    obs = obs.melt(id_vars=["SS_name", "Res_N", "Res_type"],
                   value_vars=set(obs.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    preds = preds.melt(id_vars=["Res_name", "Res_N", "Res_type", "Res_typem1"],
                   value_vars=set(preds.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    
    #Get rid of rows where the observed or predicted data is missing
    obs = obs.dropna(subset=["Shift"])
    preds = preds.dropna(subset=["Shift"])
    
    obs["ID"] = i
    preds["ID"] = i

    if type(obs_all)!=type(obs):
        obs_all = obs.copy()
        preds_all = preds.copy()
    else:
        obs_all = pd.concat([obs_all, obs], ignore_index=True)
        preds_all = pd.concat([preds_all, preds], ignore_index=True)

df = pd.merge(obs_all, preds_all, on=["ID","Res_N","Res_type", "Atom_type"], 
              how="outer", suffixes=["_obs","_pred"])

# Get rid of any residues with mismatching types, or where the type is NaN
df = df.loc[~(df["ID"]+df["Res_N"].astype(str)).isin(bad_res["ID"]+bad_res["Res_N"].astype(str)),:]
df = df.dropna(axis=0, subset=["Res_type","Res_typem1","SS_name", "Res_name"])
df["ID_Res"] = df["ID"]+ " " + df["Res_name"].astype(str)
df = df[["ID_Res", "ID","Res_name","Res_N","Res_type","Res_typem1","Atom_type","Shift_obs","Shift_pred"]]
        
# Calculate the prediction errors
df["Delta"] = df["Shift_pred"] - df["Shift_obs"]


#%% Work out obs vs pred distance matrix.
# Then calculate how many observations are closer to the prediction than the correct one

# Scaling factors
scale = (df.groupby("Atom_type").quantile(0.75) - 
         df.groupby("Atom_type").quantile(0.25))
scale = scale.loc[:, "Shift_obs"]
# Replace CB scale with CA, to avoid being skewed by Ser,Thr CB shifts.
scale["CB"] = scale["CA"]
scale["CBm1"] = scale["CAm1"]
scale = scale/scale["H"]
   
# Make a wide data frame of the obs and preds values
df_wide = df.drop(["Atom_type","Shift_obs", "Shift_pred", "Delta"], axis="columns").drop_duplicates()
df_wide.index = df_wide["ID_Res"]
df_wide = df_wide.drop("ID_Res", axis="columns")
obs_wide = pd.concat([df_wide,
                     df.pivot(index="ID_Res", 
                              columns="Atom_type", 
                              values="Shift_obs")],axis=1)
preds_wide = pd.concat([df_wide,
                     df.pivot(index="ID_Res", 
                              columns="Atom_type", 
                              values="Shift_pred")],axis=1)
  
rank_df = None

def calc_dist_matrix(x1, x2, atoms, sf, na_value=np.nan):
    """Calculate a (Euclidean) distance matrix between two sets of residues.
    
    x1 and x2: dataframes with chemical shift values in columns.
    atoms: a list or set of atoms/columns to use
    sf: a dictionary of floats, where the keys are atom names and the values 
        are the amount to scale each dimension by (eg {"H":1, "N":5} would 
        correspond to the normal values used for CSP measurements)
    na_value: if data is missing for a particular obs/pred/atom combo, what 
        value should be stored? (np.nan will give np.nan overall, 0 will 
        silently ignore the missing data)
    """
    delta2={}
    
    for atom in atoms:
        x1_atom = (x1[atom].repeat(len(x1.index)).values.
                    reshape([len(x1.index),-1]))
        x2_atom = (x2[atom].repeat(len(x2.index)).values.
                      reshape([len(x2.index),-1]).transpose())
        
        
        delta2[atom] = pd.DataFrame(((x2_atom - x1_atom)/sf[atom])**2, index=x1.index, columns=x2.index)
        # Make a note of NA positions in delta, and set them to default value 
        na_mask = np.isnan(delta2[atom])
        delta2[atom][na_mask] = na_value
    
    dist_mat = sum(delta2.values()).applymap(sqrt)
    
    return(dist_mat)
    

for id in df_wide["ID"].unique():
    print(id)
    obs_id = obs_wide.loc[(obs_wide["ID"]==id) & (obs_wide["Res_type"]!="P"),:]
    preds_id = preds_wide.loc[(preds_wide["ID"]==id) & (preds_wide["Res_type"]!="P"),:]
    obs_id.index = obs_id["Res_name"]
    preds_id.index = preds_id["Res_name"]
    
    # calculate distance between each observation and prediction
    dist_pred_HN = calc_dist_matrix(obs_id, preds_id, 
                                 atoms={"H", "N"}, sf=scale)
    dist_pred_HNCO = calc_dist_matrix(obs_id, preds_id, 
                                   atoms={"H", "N", "Cm1"}, sf=scale)
    dist_pred_CA_CO = calc_dist_matrix(obs_id, preds_id, 
                                    atoms={"H","N","C","CA","Cm1","CAm1"}, 
                                    sf=scale)
    dist_pred_most = calc_dist_matrix(obs_id, preds_id, 
                                   atoms={"H","N","C","CA","CB",
                                          "Cm1","CAm1","CBm1"}, 
                                   sf=scale)
    
    # Calculate distance between observations (penalise diagonal)
    dist_obs_HN = calc_dist_matrix(obs_id, obs_id, atoms={"H", "N"}, sf=scale) + np.diag([10]*len(obs_id.index))
    dist_obs_HNCO = calc_dist_matrix(obs_id, obs_id, atoms={"H", "N", "Cm1"}, sf=scale) + np.diag([10]*len(obs_id.index))
    dist_obs_CA_CO = calc_dist_matrix(obs_id, obs_id, atoms={"H","N","C","CA","Cm1","CAm1"}, sf=scale) + np.diag([10]*len(obs_id.index))
    dist_obs_most = calc_dist_matrix(obs_id, obs_id, atoms={"H","N","C","CA","CB","Cm1","CAm1","CBm1"}, sf=scale) + np.diag([10]*len(obs_id.index))
    
    # For each spin system, work out the rank of the correct predicted residue
    # ie. how many predictions were closer to this SS than the correct one?
    # Also store the distance to the correct prediction, the
    tmp = pd.DataFrame({"ID":id, "Res_name":dist_pred_HN.index, "N_aa":len(dist_pred_HN.index),
                        "dist_correct_most":np.diag(dist_pred_most),
                        "rank_most":np.diag(dist_pred_most.rank(axis=0)), 
                        "nearest_pred_most":dist_pred_most.min(axis=0),
                        "nearest_obs_most":dist_obs_most.min(axis=0),
                        "dist_correct_HN":np.diag(dist_pred_HN),
                        "rank_HN":np.diag(dist_pred_HN.rank(axis=0)),
                        "nearest_pred_HN":dist_pred_HN.min(axis=0),
                        "nearest_obs_HN":dist_obs_HN.min(axis=0),
                        "dist_correct_HNCO":np.diag(dist_pred_HNCO),
                        "rank_HNCO":np.diag(dist_pred_HNCO.rank(axis=0)),
                        "nearest_pred_HNCO":dist_pred_HNCO.min(axis=0),
                        "nearest_obs_HNCO":dist_obs_HNCO.min(axis=0),
                        "dist_correct_CA_CO":np.diag(dist_pred_CA_CO),
                        "rank_CA_CO":np.diag(dist_pred_CA_CO.rank(axis=0)),
                        "nearest_pred_CA_CO":dist_pred_CA_CO.min(axis=0),
                        "nearest_obs_CA_CO":dist_obs_CA_CO.min(axis=0)})
        
    if type(rank_df) != type(tmp):
        rank_df = tmp
    else:
        rank_df = pd.concat([rank_df, tmp], ignore_index=True)

#%%

# With this dataframe, we can now look at the distribution of ranks
(~rank_df.isna()).sum()     # Summarise non-na values

# Plt histograms of prediction rank
plt_base = (ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:], 
                                    mapping=aes(y="100*stat(density)")) + 
    xlim(0,50) + ylab("Percentage predictions with each rank"))
save_as_pdf_pages(filename=path/"plots/rank histograms.pdf",
                  plots=[plt_base + geom_histogram(aes(x="rank_HN"), binwidth=1),
                         plt_base + geom_histogram(aes(x="rank_HNCO"), binwidth=1),
                         plt_base + geom_histogram(aes(x="rank_CA_CO"), binwidth=1),
                         plt_base + geom_histogram(aes(x="rank_most"), binwidth=1)])

# Plot distribution of nearest neighbour vs prediction error 
plt_base = ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:])  
plt_base = plt_base + scale_colour_discrete(name="Distribution:", 
                                            labels=["Dist to correct prediction",
                                                    "Dist to nearest neighbour"])
plt_base = plt_base + theme(legend_position=(0.7,0.8)) + xlab("Scaled combined shift difference")
save_as_pdf_pages(filename=path/"plots/nearest neighbour vs prediction distributions.pdf",
                  plots=[plt_base + geom_density(aes(x="dist_correct_HN", colour=True)) + geom_density(aes(x="nearest_obs_HN", colour=False)),
                         plt_base + geom_density(aes(x="dist_correct_HNCO", colour=True)) + geom_density(aes(x="nearest_obs_HNCO", colour=False)),
                         plt_base + geom_density(aes(x="dist_correct_CA_CO", colour=True)) + geom_density(aes(x="nearest_obs_CA_CO", colour=False)),
                         plt_base + geom_density(aes(x="dist_correct_most", colour=True)) + geom_density(aes(x="nearest_obs_most", colour=False))])


    
#plt_base = ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:])
#plt_base + geom_density(aes(x="dist_correct_most", colour="ID"))
#plt_base + geom_density(aes(x="dist_correct_CA_CO"), colour="black") + geom_density(aes(x="nearest_pred_CA_CO"), colour="red")
#plt_base + geom_histogram(aes(x="dist_correct_most", y="stat(density)"), binwidth=0.1)
#
#plt_base + geom_point(aes(y="rank_HN", x="dist_correct_HNCO"), alpha=0.2)
#
#plt_base + geom_histogram(aes(x="rank_HN", y="stat(density)*100"), binwidth=1) + xlim(0,50)

#%% Plot error distribution for each atom type
plt = ggplot(data = df) + geom_density(aes(x="Delta", colour="Atom_type")) 
plt = plt + facet_wrap("Atom_type", scales="free_y")
plt = plt + xlim(-5,5) + xlab("Prediction error")
plt = plt + ggtitle("Distribution of ShiftX2 prediction errors")
plt.save(path/"plots/ShiftX2 error distribution.pdf")


plt = ggplot(data = df) + geom_point(aes(x="Shift_obs", y="Delta", colour="Res_type"))
plt = plt + facet_wrap("Atom_type", scales="free")
plt = plt + xlab("Observed shift") + ylab("Prediction error")
plt = plt + ggtitle("Prediction error is inversely correlated to chemical shift")
plt.save(path/"plots/Error vs obs shift.pdf", height=210, width=297, units="mm")

plt = ggplot(data = df) + geom_point(aes(x="Shift_pred", y="Delta", colour="Atom_type"))
plt = plt + facet_wrap("Atom_type", scales="free")
plt = plt + xlab("Predicted shift") + ylab("Prediction error")
plt.save(path/"plots/Error vs pred shift.pdf", height=210, width=297, units="mm")


#%% Make a "normalised" chemical shift that accounts for residue type and observed shift
df2 = df.copy()

df2["Res_norm"] = df2["Res_type"]
df2.loc[df2["Atom_type"].isin(["Cm1","CAm1","CBm1"]),"Res_norm"] = df2.loc[df2["Atom_type"].isin(["Cm1","CAm1","CBm1"]),"Res_typem1"]

df2["Delta_norm"] = np.nan
lm_results = pd.DataFrame(columns=["Atom_type","Res_type","Grad","Offset"])

for atom in atom_set:
    for res in df2["Res_norm"].unique():
        mask = (df2["Atom_type"]==atom) & (df2["Res_norm"]==res)
        tmp = df2.loc[mask, ["Shift_obs", "Delta"]].dropna(how="any")
        try:
            lm = linregress(tmp["Shift_obs"], tmp["Delta"])
            lm_results.loc[atom+"_"+res, :] = [atom, res, lm[0], lm[1]]
            df2.loc[mask,"Delta_norm"] = (df2.loc[mask,"Delta"] 
                                          - lm[0]*df2.loc[mask,"Shift_obs"]
                                          - lm[1])
        except:
            print("Error: ",res, atom)

plt = ggplot(data = df2.loc[df2["Res_norm"]=="N",:]) 
plt = plt + geom_point(aes(x="Shift_obs", y="Delta", colour="Res_norm"))
plt = plt + facet_wrap("Atom_type", scales="free")
#plt = plt + xlab("Observed shift") + ylab("Prediction error")
plt

#%% Look at correlated errors bertween different atom types
delta_wide = df2.pivot(index="ID_Res", columns="Atom_type", values="Delta")
delta_corr = delta_wide.corr()

# Something *really* weird is going on with the Delta values for A001 - they're 
# clustering on only a few distinct values. Maybe shiftY is involved somehow? 

ggplot(delta_wide) + geom_point(aes(x="CA", y="HA"))

#%% Not updated below this point to account for long vs wide dataframe, so won't work.


# Fit linear model between d_atom and atom_obs, and subtract to form dd_atom
atoms_m1 = pd.Series(["Cm1","CAm1","CBm1"])
lm_results = pd.DataFrame(columns=["Atom_type","Res_type","Grad","Offset"])
for atom in atom_set:
    df["dd_"+atom] = np.NaN
    
    for res in df["Res_type"].unique():
        if atom in atoms_m1:
            tmp = df.loc[df["Res_typem1"]==res, [atom+"_obs", "d_"+atom]].dropna(how="any")
        else:
            tmp = df.loc[df["Res_type"]==res, [atom+"_obs", "d_"+atom]].dropna(how="any")
        try:
            lm = linregress(tmp[atom+"_obs"], tmp["d_"+atom])
            lm_results.loc[atom+"_"+res,:] = [atom, res, lm[0], lm[1]]
            if atom in atoms_m1:
                df.loc[df["Res_typem1"]==res, "dd_"+atom] = df.loc[df["Res_typem1"]==res, "d_"+atom] - lm[0]*df.loc[df["Res_typem1"]==res, atom+"_obs"] - lm[1]
            else:
                df.loc[df["Res_type"]==res, "dd_"+atom] = df.loc[df["Res_type"]==res, "d_"+atom] - lm[0]*df.loc[df["Res_type"]==res, atom+"_obs"] - lm[1]
        except:
            print("Error: ",res, atom)
# For some reason this isn't working on the m1 atoms...

# Fit linear model between d_atom and atom_obs
#fit_results = pd.DataFrame(columns=["Res_type","Atom_type","Grad","Offset","r2","p-value","stderr"])
#for res in df["Res_type"].unique():
#    for atom in atoms_obs:
#        tmp = df.loc[df["Res_type"]==res, [atom, "d_"+atom.split("_")[0]]].dropna(how="any")
#        try:
#            lm = linregress(tmp[atom], tmp["d_"+atom.split("_")[0]])
#            fit_results.loc[res+"_"+atom,:] = [res, atom, lm[0], lm[1], lm[2]**2, lm[3], lm[4]]
#        except:
#            print("Error: ",res, atom)
#fit_results = fit_results.sort_values(by="Atom_type")
#fit_results["Grad"] = fit_results["Grad"].astype(float)
#ggplot(data=fit_results) + geom_point(aes(y="Grad", colour="Res_type", x="Atom_type")) 

#%% Calculate covariance matrixes with and without delta adjustment

d_df = df.loc[:,["d_"+ atom for atom in a.pars["atom_set"]]]
dd_df = df.loc[:,["dd_"+ atom for atom in a.pars["atom_set"]]]

d_mean = d_df.mean()
d_cov = d_df.cov()
dd_mean = dd_df.mean()
dd_cov = dd_df.cov()

# Remove "d_" or "dd_" from index and column names
d_mean.index = [s[2:] for s in d_mean.index]
d_cov.index = [s[2:] for s in d_cov.index]
d_cov.columns = [s[2:] for s in d_cov.columns]
dd_mean.index = [s[3:] for s in dd_mean.index]
dd_cov.index = [s[3:] for s in dd_cov.index]
dd_cov.columns = [s[3:] for s in dd_cov.columns]

# Save to files for use with NAPS_assigner
d_mean.to_csv("../data/d_mean.csv")
d_cov.to_csv("../data/d_cov.csv")
dd_mean.to_csv("../data/dd_mean.csv")
dd_cov.to_csv("../data/dd_cov.csv")

abs_corr_matrix = abs(dd_df.corr())

mvn = multivariate_normal(dd_df.mean(), dd_df.cov())

tmp = pd.DataFrame(mvn.rvs(size=10000), columns = dd_df.columns)
ggplot(data=tmp) + geom_point(aes(x="dd_HA", y="dd_CA"))
ggplot(data=df) + geom_point(aes(x="dd_HA", y="dd_CA"))

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
#### Exploratory analysis
plt = ggplot(data=df)
for a in ["C","CA","CB","H","HA","N","Cm1","CAm1","CBm1"]:
    plt=plt + geom_density(aes(x="d_"+a))
plt + ggtitle("Prediction error distribution (scaled)")

plt = ggplot(data=df)
for a in ["C","CA","CB","H","HA","N","Cm1","CAm1","CBm1"]:
    plt=plt + geom_density(aes(x=a+"_obs"))
plt + ggtitle("Observed shift distribution (scaled)")

ggplot(data=df) + geom_point(aes(x="CB_obs", y="d_CB", colour="Res_type")) + geom_line(aes(x='x', y='y'), data=pd.DataFrame({'x':[-7,7],'y':[0.27*7+0.01, -0.27*7+0.01]}))

#ggplot(data=fit_results) + geom_point(aes(y="Grad", colour="Res_type", x="Atom_type")) 

ggplot(data=df) + geom_point(aes(x="N_obs", y="dd_N", colour="N_pred"))

atom = "N"
plt = ggplot(data=df) + geom_density(aes(x=atom+"_obs"))
lo = df[atom+"_obs"].min()
hi = df[atom+"_obs"].max()
step = (hi-lo)/6
plt = ggplot(aes(x="dd_"+atom)) + geom_density(data=df.loc[df[atom+"_obs"]<lo+step,:], colour="red") 
plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+step,lo+2*step),:], colour="orange")
plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+2*step,lo+3*step),:], colour="yellow")
plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+3*step,lo+4*step),:], colour="green")
plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+4*step,lo+5*step),:], colour="blue")
plt = plt + geom_density(data=df.loc[df[atom+"_obs"]>lo+5*step,:], colour="magenta")
plt

ggplot(data=df) + geom_density(aes(x="CA_obs", colour="Res_type"))
ggplot(data=df) +geom_point(aes(x="CA_obs", y="d_CA", colour="Res_type"), alpha=0.2) + facet_wrap("Res_type")

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



ggplot(data=df) + geom_point(aes(x="CA_obs", y="d_CA")) + geom_smooth(aes(x="CA_obs", y="d_CA"), colour="red", method="mavg")


ggplot(data=rank_df) + geom_bar(aes(x="rank_all"), stat=stat_bin(binwidth=1))
