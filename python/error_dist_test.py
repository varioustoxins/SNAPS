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

path = Path("..")
#path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/Alex/GitHub/NAPS/")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")
#%% Prepare to import all proteins in the testset

testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, 
                           names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/"data/testset/simplified_BMRB"/file 
                      for file in testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/"data/testset/shiftx2_results"/file 
                      for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df.index = testset_df["ID"]

importer = NAPS_importer()
assigner = NAPS_assigner()

# Work out which residue types don't match between obs and preds
seq_df = None
for i in testset_df["ID"]:
    print(i)
    obs = importer.import_testset_shifts(testset_df.loc[i, "obs_file"])
    # Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    
#    # Make columns for the i-1 predicted shifts of C, CA and CB
#    obs.index = obs["Res_N"]
#    obs.index.name = None
#    obs_m1 = obs[["Res_type","Res_N"]].copy()
#    obs_m1.index = obs_m1.index+1
#    obs_m1.columns = obs_m1.columns + "_m1"
#    obs = pd.merge(obs, obs_m1, how="left", left_index=True, right_index=True)
#    obs.index = obs["SS_name"]
#    obs.index.name = None
    
    preds = assigner.import_pred_shifts(testset_df.loc[i, "preds_file"], 
                                        filetype="shiftx2")
    #Change Bs to Cs
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type_m1"]=="B", "Res_type_m1"] = "C"
       
    df = pd.merge(obs, preds, on="Res_N", how="outer", suffixes=["_obs","_pred"])
    df["ID"] = i
    if type(seq_df)!=type(df):
        seq_df = df.copy()
    else:
        seq_df = pd.concat([seq_df, df], ignore_index=True)

seq_df = seq_df[["ID","Res_N","Res_type_obs","Res_type_pred",
                 "Res_type_m1_obs","Res_type_m1_pred"]]
seq_df = seq_df.fillna("NA")
tmp1 = seq_df.loc[seq_df["Res_type_obs"]!=seq_df["Res_type_pred"],:]
tmp1 = tmp1.loc[~tmp1["Res_type_obs"].isin(["NA"]) & ~tmp1["Res_type_pred"].isin(["NA"]),:]
tmp2 = seq_df.loc[seq_df["Res_type_m1_obs"]!=seq_df["Res_type_m1_pred"],:]
tmp2 = tmp2.loc[~tmp2["Res_type_m1_obs"].isin(["NA"]) & ~tmp2["Res_type_m1_pred"].isin(["NA"]),:]

bad_res = pd.concat([tmp1, tmp2])
bad_res = bad_res.drop_duplicates()

#%% Import the actual data

importer = NAPS_importer()
assigner = NAPS_assigner()

atom_set = {"H","N","HA","C","CA","CB","C_m1","CA_m1","CB_m1"}

obs_all = None
for i in testset_df["ID"]:
#for i in ["A063"]:
    print(i)
    obs = importer.import_testset_shifts(testset_df.loc[i, "obs_file"], 
                                         remove_Pro=False)
    preds = assigner.import_pred_shifts(testset_df.loc[i, "preds_file"], 
                                        filetype="shiftx2")
    
    #Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type"]=="B", "Res_name"] = (
            preds.loc[preds["Res_type"]=="B", "Res_name"].str.replace("B","C"))
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type_m1"]=="B", "Res_type_m1"] = "C"
    

    #Convert wide to long
    obs = obs.melt(id_vars=["SS_name", "Res_N", "Res_type"],
                   value_vars=set(obs.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    preds = preds.melt(id_vars=["Res_name", "Res_N", "Res_type", "Res_type_m1"],
                   value_vars=set(preds.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    
    #Get rid of rows where the observed or predicted data is missing
    #obs = obs.dropna(subset=["Shift"])
    #preds = preds.dropna(subset=["Shift"])
    
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
df = df.dropna(axis=0, subset=["Res_type","Res_type_m1","SS_name", "Res_name"])
df["ID_Res"] = df["ID"]+ " " + df["Res_name"].astype(str)
df = df[["ID_Res", "ID","Res_name","Res_N","Res_type","Res_type_m1",
         "Atom_type","Shift_obs","Shift_pred"]]
        
# Calculate the prediction errors
df["Delta"] = df["Shift_pred"] - df["Shift_obs"]

df_keep_na = df.copy()
df = df.dropna(axis=0, subset=["Delta"])

# Save the data
df.to_csv(path/"output/prediction_accuracy.txt", sep="\t", float_format="%.3f")

#%% Work out obs vs pred distance matrix.
# Then calculate how many observations are closer to the prediction than the correct one

# Scaling factors
df_q = (df.groupby("Atom_type").quantile(
        [0,0.025,0.1,0.25,0.5,0.75,0.9,0.975,1])["Shift_obs"].unstack())

# Recalculate for CB, removing S and T
df_q.loc["CB",:] = df[(df["Atom_type"]=="CB") & (~df["Res_type"].isin(["S","T"]))].quantile(
        [0,0.025,0.1,0.25,0.5,0.75,0.9,0.975,1])["Shift_obs"]
df_q.loc["CB_m1",:] = df[(df["Atom_type"]=="CB_m1") & (~df["Res_type_m1"].isin(["S","T"]))].quantile(
        [0,0.025,0.1,0.25,0.5,0.75,0.9,0.975,1])["Shift_obs"]


df_range = df_q[1] - df_q[0]
df_95 = df_q[0.975] - df_q[0.025]
df_80 = df_q[0.9] - df_q[0.1]
df_IQR = df_q[0.75] - df_q[0.25]

scale = (df.groupby("Atom_type").quantile(0.75) - 
         df.groupby("Atom_type").quantile(0.25))
scale = scale.loc[:, "Shift_obs"]
# Replace CB scale with CA, to avoid being skewed by Ser,Thr CB shifts.
scale["CB"] = scale["CA"]
scale["CB_m1"] = scale["CA_m1"]
scale = scale/scale["H"]
   
# Make a wide data frame of the obs and preds values
df_wide = df.drop(["Atom_type","Shift_obs", "Shift_pred", "Delta"], 
                  axis="columns").drop_duplicates()
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
        
        
        delta2[atom] = pd.DataFrame(((x2_atom - x1_atom)/sf[atom])**2, 
                                    index=x1.index, columns=x2.index)
        # Make a note of NA positions in delta, and set them to default value 
        na_mask = np.isnan(delta2[atom])
        if atom in ("CB","CB_m1", "HA"):    
            # Glycines don't have CB and HA, so treat these atoms specially
            delta2[atom][na_mask] = 0
        else:
            delta2[atom][na_mask] = na_value
    
    dist_matrix = sum(delta2.values()).applymap(sqrt)
    
    return(dist_matrix)
    
def penalise_aa_mismatch(x1, x2, dist_matrix, penalty=1000, 
                         aa_groups=None, aa_groups_m1=None):
    """
    Penalise cells in a distance matrix where amino acid types don't match.
    
    aa_groups: list of groups of amino acids (eg. 
        ["A","S","T","G","CDEFHIKLMNPQRVWY"]. Any pairs between x1 and x2 which 
        aren't part of the same group will be penalised. If None, no penalties 
        apply.
    aa_groups_m1: As above, but bsed on i-1 residue
    """
    if aa_groups!=None:
        penalty_matrix = pd.DataFrame(penalty, index=dist_matrix.index, 
                                      columns=dist_matrix.columns)
        for aa_group in aa_groups:
            penalty_matrix.loc[x1["Res_type"].isin(list(aa_group)), 
                            x2["Res_type"].isin(list(aa_group))] = 0          
        dist_matrix = dist_matrix + penalty_matrix
        
    if aa_groups_m1!=None:
        penalty_matrix = pd.DataFrame(penalty, index=dist_matrix.index, 
                                      columns=dist_matrix.columns)
        for aa_group in aa_groups_m1:
            penalty_matrix.loc[x1["Res_type_m1"].isin(list(aa_group)), 
                            x2["Res_type_m1"].isin(list(aa_group))] = 0          
        dist_matrix = dist_matrix + penalty_matrix
    
    return(dist_matrix)

for id in df_wide["ID"].unique():
    print(id)
    obs_id = obs_wide.loc[(obs_wide["ID"]==id) & (obs_wide["Res_type"]!="P"),:]
    preds_id = preds_wide.loc[(preds_wide["ID"]==id) & (preds_wide["Res_type"]!="P"),:]
    obs_id.index = obs_id["Res_name"]
    preds_id.index = preds_id["Res_name"]
    
    # calculate distance between each observation and prediction
    dist_pred_HN = calc_dist_matrix(obs_id, preds_id, sf=scale, 
                                 atoms={"H", "N"})
    dist_pred_HNCO = calc_dist_matrix(obs_id, preds_id, sf=scale, 
                                   atoms={"H", "N", "C_m1"})
    dist_pred_CA_CO = calc_dist_matrix(obs_id, preds_id, sf=scale,
                                    atoms={"H","N","C","CA","C_m1","CA_m1"})
    dist_pred_most = calc_dist_matrix(obs_id, preds_id, sf=scale,
                                   atoms={"H","N","C","CA","CB",
                                          "C_m1","CA_m1","CB_m1"})
    dist_pred_HN_HADAMAC = penalise_aa_mismatch(obs_id, preds_id, dist_pred_HN,
                                                aa_groups_m1=["S","T","G","AVI",
                                                      "DN","FHYWC","REKPQML"])
    dist_pred_HNCO_HADAMAC = penalise_aa_mismatch(obs_id, preds_id, dist_pred_HNCO,
                                            aa_groups_m1=["S","T","G","AVI",
                                                    "DN","FHYWC","REKPQML"])
    
    # Calculate distance between observations (penalise diagonal)
    dist_obs_HN = (calc_dist_matrix(obs_id, obs_id, sf=scale,
                                    atoms={"H", "N"}) + 
                   np.diag([1000]*len(obs_id.index)))
    dist_obs_HNCO = (calc_dist_matrix(obs_id, obs_id, sf=scale,
                                      atoms={"H", "N", "C_m1"}) + 
                     np.diag([1000]*len(obs_id.index)))
    dist_obs_CA_CO = (calc_dist_matrix(obs_id, obs_id, sf=scale,
                                       atoms={"H","N","C","CA","C_m1","CA_m1"}) + 
                      np.diag([1000]*len(obs_id.index)))
    dist_obs_most = (calc_dist_matrix(obs_id, obs_id, sf=scale,
                                      atoms={"H","N","C","CA","CB","C_m1",
                                             "CA_m1","CB_m1"}) + 
                     np.diag([1000]*len(obs_id.index)))
    dist_obs_HN_HADAMAC = penalise_aa_mismatch(obs_id, obs_id, dist_obs_HN, 
                                            aa_groups_m1=["S","T","G","AVI",
                                                    "DN","FHYWC","REKPQML"])
    dist_obs_HNCO_HADAMAC = penalise_aa_mismatch(obs_id, obs_id, dist_obs_HNCO,
                                            aa_groups_m1=["S","T","G","AVI",
                                                    "DN","FHYWC","REKPQML"])
    
    # For each spin system, work out the rank of the correct predicted residue
    # ie. how many predictions were closer to this SS than the correct one?
    # Also store the distance to the correct prediction, the
    tmp = pd.DataFrame({"ID":id, "Res_name":dist_pred_HN.index, 
                        "N_aa":len(dist_pred_HN.index),
                        "dist_correct_most":np.diag(dist_pred_most),
                        "rank_most":np.diag(dist_pred_most.rank(axis=0)), 
                        "nearest_wrong_pred_most":(dist_pred_most+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_most":dist_obs_most.min(axis=0),
                        "dist_correct_HN":np.diag(dist_pred_HN),
                        "rank_HN":np.diag(dist_pred_HN.rank(axis=0)),
                        "nearest_wrong_pred_HN":(dist_pred_HN+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_HN":dist_obs_HN.min(axis=0),
                        "dist_correct_HNCO":np.diag(dist_pred_HNCO),
                        "rank_HNCO":np.diag(dist_pred_HNCO.rank(axis=0)),
                        "nearest_wrong_pred_HNCO":(dist_pred_HNCO+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_HNCO":dist_obs_HNCO.min(axis=0),
                        "dist_correct_CA_CO":np.diag(dist_pred_CA_CO),
                        "rank_CA_CO":np.diag(dist_pred_CA_CO.rank(axis=0)),
                        "nearest_wrong_pred_CA_CO":(dist_pred_CA_CO+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_CA_CO":dist_obs_CA_CO.min(axis=0),
                        "dist_correct_HN_HADAMAC":np.diag(dist_pred_HN_HADAMAC),
                        "rank_HN_HADAMAC":np.diag(dist_pred_HN_HADAMAC.rank(axis=0)),
                        "nearest_wrong_pred_HN_HADAMAC":(dist_pred_HN_HADAMAC+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_HN_HADAMAC":dist_obs_HN_HADAMAC.min(axis=0),
                        "dist_correct_HNCO_HADAMAC":np.diag(dist_pred_HNCO_HADAMAC),
                        "rank_HNCO_HADAMAC":np.diag(dist_pred_HNCO_HADAMAC.rank(axis=0)),
                        "nearest_wrong_pred_HNCO_HADAMAC":(dist_pred_HNCO_HADAMAC+np.diag([1000]*len(obs_id.index))).min(axis=0),
                        "nearest_obs_HNCO_HADAMAC":dist_obs_HNCO_HADAMAC.min(axis=0)})
        
    if type(rank_df) != type(tmp):
        rank_df = tmp
    else:
        rank_df = pd.concat([rank_df, tmp], ignore_index=True)

#%% Look at distribution of nearest neighbours and ranks

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

# Plot distributions of prediction error and nearest incorrect predicition
def nearest_neighbour_dist_plots():
    plt_base = ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:])  
    plt_base = plt_base + scale_colour_brewer(name="Distribution:", 
                                              type="qual", palette=6,
                                                labels=["Dist to correct prediction",
                                                        "Dist to nearest wrong prediction"])
    plt_base = plt_base  + theme_bw()+ theme(legend_position="bottom")
    plt_base = plt_base + xlab("Scaled combined shift difference") + xlim(0,3)
    
    for x in ["HN","HNCO","HN_HADAMAC","HNCO_HADAMAC","CA_CO", "most"]:
        plt = plt_base
        plt = plt + geom_density(aes(x="dist_correct_"+x, colour="str(1)")) 
        plt = plt + geom_density(aes(x="nearest_wrong_pred_"+x, colour="str(2)"))
        #plt = plt + geom_density(aes(x="nearest_obs_"+x, colour="str(3)"))
        plt = plt + ggtitle("Distance distributions for: "+x)
        yield(plt)
    

save_as_pdf_pages(filename=path/"plots/nearest prediction distribution.pdf",
                  plots=nearest_neighbour_dist_plots())

#### Make figures for NMR-DG 2019 poster
# Ranks
tmp = rank_df.loc[~rank_df["dist_correct_most"].isna(),:]
tmp.loc[tmp["rank_most"]>=20, "rank_most"] = 20
tmp.loc[tmp["rank_HN"]>=20, "rank_HN"] = 20
plt_base = ggplot(data=tmp, mapping=aes(y="100*stat(density)")) 
plt_base += xlim(0,21) 
plt_base += ylim(0,100)
plt_base += ylab("Percentage of predictions with each rank")
plt_base += xlab("Rank of the correct prediction")
plt_base += theme_bw()

plt1 = plt_base + geom_histogram(aes(x="rank_HN"), breaks=range(21))
plt1 += ggtitle("Rank of correct predictions: HN")
plt1.save(path/"plots/Poster rank dist HN.pdf", height=100, width=100, units="mm")

plt2 = plt_base + geom_histogram(aes(x="rank_most"), breaks=range(21))
plt2 += ggtitle("Rank of correct predictions: most")
plt2.save(path/"plots/Poster rank dist most.pdf", height=100, width=100, units="mm")
 
# Distributions
plt_base = ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:])  
plt_base = plt_base + scale_fill_brewer(name="Distribution:", 
                                          type="qual", palette=6,
                                            labels=["Dist to correct prediction",
                                                    "Dist to nearest wrong prediction"])
plt_base = plt_base  + theme_bw()+ theme(legend_position="bottom") +ylim(0,7.5)
plt_base = plt_base + xlab("Scaled combined shift difference") + xlim(0,2)

plt1 = plt_base
plt1 = plt1 + geom_density(aes(x="dist_correct_HN", fill="str(1)"), alpha=0.7) 
plt1 = plt1 + geom_density(aes(x="nearest_wrong_pred_HN", fill="str(2)"), alpha=0.7)
plt1 = plt1 + ggtitle("Distance distributions for: HN")
plt1.save(path/"plots/Poster prediction dist HN.pdf", height=100, width=100, units="mm")

plt2 = plt_base
plt2 = plt2 + geom_density(aes(x="dist_correct_most", fill="str(1)"), alpha=0.7) 
plt2 = plt2 + geom_density(aes(x="nearest_wrong_pred_most", fill="str(2)"), alpha=0.7)
plt2 = plt2 + ggtitle("Distance distributions for: most")
plt2.save(path/"plots/Poster prediction dist most.pdf", height=100, width=100, units="mm")

#plt_base = ggplot(data=rank_df.loc[~rank_df["dist_correct_most"].isna(),:])
#plt_base + geom_density(aes(x="dist_correct_most", colour="ID"))
#plt_base + geom_density(aes(x="dist_correct_CA_CO"), colour="black") + geom_density(aes(x="nearest_pred_CA_CO"), colour="red")
#plt_base + geom_histogram(aes(x="dist_correct_most", y="stat(density)"), binwidth=0.1)
#
#plt_base + geom_point(aes(y="rank_HN", x="dist_correct_HNCO"), alpha=0.2)
#
#plt_base + geom_histogram(aes(x="rank_HN", y="stat(density)*100"), binwidth=1) + xlim(0,50)

#%% Plot error distribution for each atom type
# Calculate the overall standard deviation for each residue
delta_stdev = df.groupby("Atom_type").std()["Delta"]

# Make a column associating the i-1 residue type with i-1 shifts.
df["Res_atom"] = df["Res_type"]
df.loc[df["Atom_type"].isin(["C_m1","CA_m1","CB_m1"]),"Res_atom"] = (
         df.loc[df["Atom_type"].isin(["C_m1","CA_m1","CB_m1"]),"Res_type_m1"])

    
plt = ggplot(data = df[~df["Atom_type"].isin(["C_m1","CA_m1","CB_m1"])]) 
plt = plt + geom_density(aes(x="Delta", fill="Atom_type")) 
plt = plt + facet_wrap("Atom_type", scales="free_y")
plt = plt + xlim(-5,5) + xlab("Prediction error (ppm)")
plt = plt + ggtitle("Distribution of ShiftX2 prediction errors")
plt = plt + theme_bw()
plt = plt + scale_fill_brewer(name="Atom", type="qual", palette=6)
plt.save(path/"plots/ShiftX2 error distribution.pdf")
plt.save(path/"plots/Poster shiftX2 error distribution.pdf", height=100, width=150, units="mm")

plt = ggplot(data = df) 
plt = plt + geom_point(aes(x="Shift_obs", y="Delta", colour="Res_atom"))
plt = plt + scale_x_reverse()
plt = plt + facet_wrap("Atom_type", scales="free")
plt = plt + xlab("Observed shift") + ylab("Prediction error")
plt = plt + ggtitle("Prediction error is inversely correlated to chemical shift")
plt.save(path/"plots/Error vs obs shift.pdf", height=210, width=297, units="mm")

plt = ggplot(data = df)
plt = plt + geom_point(aes(x="Shift_pred", y="Delta", colour="Res_atom"))
plt = plt + scale_x_reverse()
plt = plt + facet_wrap("Atom_type", scales="free")
plt = plt + xlab("Predicted shift") + ylab("Prediction error")
plt.save(path/"plots/Error vs pred shift.pdf", height=210, width=297, units="mm")

plt = ggplot(data=df)
plt = plt + geom_density(aes(x="Delta", colour="ID"))
plt = plt + facet_wrap("Atom_type", scales="free")
plt.save(path/"plots/Error distribution per ID.pdf", height=210, width=297, units="mm")


def error_boxplots():
    for a in atom_set:
        plt = ggplot(df[df["Atom_type"]==a])
        plt = plt + geom_boxplot(aes(y="Delta", x="ID"))
        plt = plt + theme(axis_text_x = element_text(angle=90))
        plt = plt + ggtitle("Error distribution for atom "+a)
        yield(plt)
save_as_pdf_pages(filename=path/"plots/Error distribution per atom.pdf",
                  plots=error_boxplots())

# Plot the error standard deviation for each ID, in case some are less accurate
tmp = df.groupby(["ID", "Atom_type"])["Delta"].std().reset_index()
plt = ggplot(tmp) + geom_point(aes(x="ID", y="Delta", colour="Atom_type"))
plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=0.5))
plt.save(path/"plots/Error stdev per ID.pdf", height=210, width=297, units="mm")

# Plot how much data is missing
sum(df_keep_na["Shift_obs"].isna())
df_keep_na["Missing_obs"] = df_keep_na["Shift_obs"].isna()
tmp = df_keep_na.groupby(["ID","Atom_type"])["Missing_obs"].sum().reset_index()
plt = ggplot(tmp) + geom_point(aes(x="ID", y="Missing_obs", colour="Atom_type"))
plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=0.5))
plt.save(path/"plots/Missing obs per ID.pdf", height=210, width=297, units="mm")


#%% Make a "corrected" predicted shift that accounts for residue type and observed shift

df2 = df.copy()

df2["Res_atom"] = df2["Res_type"]
df2.loc[df2["Atom_type"].isin(["C_m1","CA_m1","CB_m1"]),"Res_atom"] = (
         df2.loc[df2["Atom_type"].isin(["C_m1","CA_m1","CB_m1"]),"Res_type_m1"])

df2["Delta_cor"] = np.nan
lm_results = pd.DataFrame(columns=["Atom_type","Res_type","Grad","Offset"])

for atom in atom_set:
    for res in df2["Res_atom"].unique():
        mask = (df2["Atom_type"]==atom) & (df2["Res_atom"]==res)
        tmp = df2.loc[mask, ["Shift_obs", "Delta"]].dropna(how="any")
        try:
            lm = linregress(tmp["Shift_obs"], tmp["Delta"])
            lm_results.loc[atom+"_"+res, :] = [atom, res, lm[0], lm[1]]
            df2.loc[mask,"Shift_pred_cor"] = (df2.loc[mask,"Shift_pred"] 
                                              - lm[0]*df2.loc[mask,"Shift_obs"]
                                              - lm[1])
            df2.loc[mask,"Delta_cor"] = (df2.loc[mask,"Shift_pred_cor"] - 
                                         df2.loc[mask,"Shift_obs"])
        except:
            print("Error: ",res, atom)

lm_results.to_csv("../config/lin_model_shiftx2.csv")

df2["Delta2"] = abs(df2["Delta_cor"]) - abs(df2["Delta"])

#%% Make some graphs showing that the corrected predictions are better

plt = ggplot(data = df2) 
plt = plt + geom_density(aes(x="Delta_cor", 
                             linetype=["Corrected"], 
                             colour="Atom_type"))
plt = plt + geom_density(aes(x="Delta", 
                             linetype=["Original"], 
                             colour="Atom_type")) 
plt = plt + facet_wrap("Atom_type", scales="free_y")
plt = plt + xlim(-5,5) + xlab("Prediction error")
plt = plt + ggtitle("Corrected prediction error has a slightly tighter "+
                    "distribution than the original predictions")
plt = plt + scale_linetype_discrete(name="Prediction type")
plt.save(path/"plots/Corrected error distribution.pdf", 
         height=210, width=297, units="mm")

plt1 = ggplot(data = df2)
plt1 = plt1 + geom_point(aes(x="Shift_obs", y="Delta_cor", colour="Res_atom"))
plt1 = plt1 + scale_x_reverse()
plt1 = plt1 + facet_wrap("Atom_type", scales="free")
plt1 = plt1 + xlab("Observed shift") + ylab("Corrected prediction error")
plt1 = plt1 + ggtitle("Corrected prediction error is independent of observed shift...")

plt2 = ggplot(data = df2)
plt2 = plt2 + geom_point(aes(x="Shift_pred", y="Delta_cor", colour="Res_atom"))
plt2 = plt2 + scale_x_reverse()
plt2 = plt2 + facet_wrap("Atom_type", scales="free")
plt2 = plt2 + xlab("Predicted shift") + ylab("Corrected prediction error")
plt2 = plt2 + ggtitle("...but is dependent on predicted shift. But I don't think this is a problem.")

save_as_pdf_pages(filename=path/"plots/Corrected delta vs shift.pdf",
                  plots=[plt1, plt2])

plt = ggplot(data=df2) + geom_density(aes(x="Delta2")) 
plt = plt + facet_wrap("Atom_type", scales="free")
plt = plt + ggtitle("Improvement in corrected Delta over original Delta")
plt = plt +  xlab("abs(Delta_cor) - abs(Delta)")
plt = plt + geom_vline(xintercept=0, colour="red")
plt.save(filename=path/"plots/Improvement in corrected Delta.pdf", 
         height=210, width=297, units="mm")

plt = ggplot(data=df2[df2["Res_atom"]=="L"]) 
plt = plt + geom_point(aes(x="Delta", y="Delta2", colour="Res_atom"))
plt = plt + facet_wrap("Atom_type", scales="free")
plt.save(filename=path/"plots/test.pdf", 
         height=210, width=297, units="mm")

plt = ggplot(data=df2[df2["Res_atom"]=="L"]) 
plt = plt + geom_point(aes(x="Shift_obs", y="Delta2"))#, colour="Res_atom"))
plt = plt + facet_wrap("Atom_type", scales="free")
plt.save(filename=path/"plots/test2.pdf", 
         height=210, width=297, units="mm")


def dist_comparison_plots():
    for res in df2["Res_atom"].unique():
        plt = ggplot(df2.loc[df2["Res_atom"]==res,:]) 
        plt = plt + geom_density(aes(x="Shift_obs", colour=["Obs"]))
        plt = plt + geom_density(aes(x="Shift_pred", colour=["Pred"]))
        plt = plt + geom_density(aes(x="Shift_pred-Delta+Delta_cor", colour=["Pred_corr"]))
        plt = plt + scale_x_reverse()
        plt = plt + facet_wrap("Atom_type", scales="free")
        plt = plt + ggtitle("Shift distribution for residue "+res)
        yield(plt)

save_as_pdf_pages(filename=path/"plots/Obs pred shift distributions.pdf", 
                  plots=dist_comparison_plots())

#%% Look at correlated errors between different atom types
delta_wide = df2.pivot(index="ID_Res", columns="Atom_type", values="Delta")
delta_corr = delta_wide.corr()

delta_cor_wide = df2.pivot(index="ID_Res", columns="Atom_type", values="Delta_cor")
delta_cor_corr = delta_cor_wide.corr()

# Save to files for use with NAPS_assigner
delta_wide.mean().to_csv("../config/d_mean.csv")
delta_wide.cov().to_csv("../config/d_cov.csv")
delta_cor_wide.mean().to_csv("../config/dd_mean.csv")
delta_cor_wide.cov().to_csv("../config/dd_cov.csv")

obs_wide = df2.pivot(index="ID_Res", columns="Atom_type", values="Shift_obs")
obs_corr = obs_wide.corr()

# Something *really* weird is going on with the Delta values for A001 - they're 
# clustering on only a few distinct values. Maybe shiftY is involved somehow? 

ggplot(delta_wide) + geom_point(aes(x="CA", y="HA"))

#%% Find all ID's which contain particular shifts
tmp=df.groupby(["ID","Atom_type"]).size().unstack(fill_value=0)
#tmp.index[tmp["HA"]==0]
#tmp.index[tmp["C"]==0]
#tmp.index[tmp["C_m1"]==0]
#tmp.index[tmp["CB"]==0]
#tmp.index[tmp["CB_m1"]==0]
# ID's which have at least some observations for all three carbons
tmp.index[(tmp["CA"]>0) & (tmp["CB"]>0) & (tmp["C"]>0)]

#%% Check for regions which lack either obs or preds
tmp = obs_all.groupby("ID")["Res_N"]
tmp2 = pd.DataFrame({"Obs_min":tmp.min(), "Obs_max":tmp.max()})
tmp = preds_all.groupby("ID")["Res_N"]
tmp2["Preds_min"] = tmp.min()
tmp2["Preds_max"] = tmp.max()
tmp2["Common_length"] = tmp2[["Obs_max","Preds_max"]].min(axis=1) - tmp2[["Obs_min","Preds_min"]].max(axis=1)
tmp2["Overhang"] = (tmp2[["Obs_max","Preds_max"]].max(axis=1) - tmp2[["Obs_min","Preds_min"]].min(axis=1)) - tmp2["Common_length"]

#%% Not updated below this point to account for long vs wide dataframe, so won't work.




#%% Calculate covariance matrixes with and without delta adjustment

#abs_corr_matrix = abs(dd_df.corr())
#
#mvn = multivariate_normal(dd_df.mean(), dd_df.cov())
#
#tmp = pd.DataFrame(mvn.rvs(size=10000), columns = dd_df.columns)
#ggplot(data=tmp) + geom_point(aes(x="dd_HA", y="dd_CA"))
#ggplot(data=df) + geom_point(aes(x="dd_HA", y="dd_CA"))

#%%
#### Try out multivariate normal distribution
#df2 = df[["ID","SS_name","Res_N","Res_type","Res_type_m1"]+[s+"_obs" for s in a.pars["atom_set"]]+["d_"+s for s in a.pars["atom_set"]]]
#
#b = abs(df2.corr())
#b[b<0.1]=-1
#
#df2_mean = df2.mean()
#df2_mean = df2_mean.drop("Res_N")
#df2_cov = df2.cov()
#df2_cov = df2_cov.drop(index="Res_N", columns="Res_N")
##df2_cov[abs(df2_cov)<0.1] = 0 # Can't set small covariances to zero because it makes matrix non-invertible, which causes error in multivariate_normal()
#
#
#mvn_full = multivariate_normal(df2_mean, df2_cov)
#mvn_obs = multivariate_normal(df2_mean.iloc[0:9], df2_cov.iloc[0:9, 0:9])   # The expected distribution of the observations
#mvn_d = multivariate_normal(df2_mean.iloc[9:18], df2_cov.iloc[9:18, 9:18])  # This is the marginal distribution of the errors over the observations
#
## Test random variables drawn from distribution
#tmp = pd.DataFrame(mvn_full.rvs(size=10000), columns = df2_cov.columns)
#ggplot(data=tmp) + geom_point(aes(x="d_CB_m1", y="CB_m1_obs"))    # Distribution is more compressed than real one. I guess this is because the real distribution has longer tails than the normal
#c = tmp.cov()
#d = df2_cov - c     # The covariance of the results matches the original data fairly well at least
#
## How to add in conditional probability?
#obs1 = [0]*18
#mvn_full.pdf(obs1)
#
#mvn_d.pdf(obs1[9:18])
#
#from math import exp
#e = exp(mvn_full.logpdf(obs1) - mvn_obs.logpdf(obs1[0:9]))


#%%
#### Exploratory analysis
#plt = ggplot(data=df)
#for a in ["C","CA","CB","H","HA","N","C_m1","CA_m1","CB_m1"]:
#    plt=plt + geom_density(aes(x="d_"+a))
#plt + ggtitle("Prediction error distribution (scaled)")
#
#plt = ggplot(data=df)
#for a in ["C","CA","CB","H","HA","N","C_m1","CA_m1","CB_m1"]:
#    plt=plt + geom_density(aes(x=a+"_obs"))
#plt + ggtitle("Observed shift distribution (scaled)")
#
#ggplot(data=df) + geom_point(aes(x="CB_obs", y="d_CB", colour="Res_type")) + geom_line(aes(x='x', y='y'), data=pd.DataFrame({'x':[-7,7],'y':[0.27*7+0.01, -0.27*7+0.01]}))
#
##ggplot(data=fit_results) + geom_point(aes(y="Grad", colour="Res_type", x="Atom_type")) 
#
#ggplot(data=df) + geom_point(aes(x="N_obs", y="dd_N", colour="N_pred"))
#
#atom = "N"
#plt = ggplot(data=df) + geom_density(aes(x=atom+"_obs"))
#lo = df[atom+"_obs"].min()
#hi = df[atom+"_obs"].max()
#step = (hi-lo)/6
#plt = ggplot(aes(x="dd_"+atom)) + geom_density(data=df.loc[df[atom+"_obs"]<lo+step,:], colour="red") 
#plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+step,lo+2*step),:], colour="orange")
#plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+2*step,lo+3*step),:], colour="yellow")
#plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+3*step,lo+4*step),:], colour="green")
#plt = plt + geom_density(data=df.loc[df[atom+"_obs"].between(lo+4*step,lo+5*step),:], colour="blue")
#plt = plt + geom_density(data=df.loc[df[atom+"_obs"]>lo+5*step,:], colour="magenta")
#plt
#
#ggplot(data=df) + geom_density(aes(x="CA_obs", colour="Res_type"))
#ggplot(data=df) +geom_point(aes(x="CA_obs", y="d_CA", colour="Res_type"), alpha=0.2) + facet_wrap("Res_type")
#
#ggplot(data=df2) + geom_point(aes(x="d_CB_m1", y="CB_m1_obs", colour="Res_type_m1"))# + facet_wrap("Res_type")
##ggplot(data=df[df["Res_type"]!="C"])
#
#tmp = df2.loc[df2["CB_m1_obs"]>7.5,:]
#
#
#leu = df.loc[df["Res_type"]=="L", ["d_"+atom for atom in a.pars["atom_set"]]]
#
##ggplot(data=leu) + geom_point(aes(x="CB_obs", y="CB_pred"))
#
#leu_mean = leu.mean(axis=0)
#leu_cov = leu.cov()
#leu_corr = leu.corr()
#
#delta = df.loc[:, ["ID","Res_name","Res_type", "Res_type_m1"]+["d_"+atom for atom in a.pars["atom_set"]]]
#d_mean = delta.mean(axis=0)
#d_cov = delta.cov()
#d_corr = delta.corr()
#
#ggplot(data=delta) + geom_point(aes(x="d_CA_m1", y="d_CB_m1", colour="Res_type"))
#
#
#
#tmp = delta.loc[abs(delta["d_CB_m1"])>5,:]
#tmp2 = df.loc[(df["ID"]+df["Res_name"]).isin(tmp["ID"]+tmp["Res_name"]),:]
#
#
#
#ggplot(data=df) + geom_point(aes(x="CA_obs", y="d_CA")) + geom_smooth(aes(x="CA_obs", y="d_CA"), colour="red", method="mavg")
#
#
#ggplot(data=rank_df) + geom_bar(aes(x="rank_all"), stat=stat_bin(binwidth=1))
