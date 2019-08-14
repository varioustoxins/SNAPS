#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Defines functions useful for analysing NAPs data, or used for testing purposes.

@author: aph516
"""

import pandas as pd
import numpy as np

def import_testset_metadata(path):
    """Import metadata on the test datasets
    
    path: top level NAPS directory
    """
    testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, 
                               names=["ID","PDB","BMRB","Resolution","Length"])
    testset_df["obs_file"] = [path/x for x in "data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"]
    testset_df["preds_file"] = [path/x for x in "data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
    testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
    testset_df.index = testset_df["ID"]
    
    return(testset_df)
    
def check_assignment_accuracy(data_dir, testset_df, ID_list, 
                              ranks=[1], prefix=""):
    """Function to check assignment accuracy
    
    Returns a list with two items: the complete set of assignments for all 
    proteins, and a summary of the accuracy for each protein
    """
    # Nb. make sure data_dir ends with a forward slash

    assigns = None
    for i in ID_list:
        tmp_all = pd.read_csv(
                data_dir/(prefix+testset_df.loc[i, "out_name"]+".txt"), 
                sep="\t", index_col=False)
        tmp_all["ID"] = i
        
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp_all["Res_type"]=="B"
        tmp_all.loc[mask,"Res_type"] = "C"
        tmp_all.loc[mask, "Res_name"] = (tmp_all.loc[mask,"Res_name"].
                   str.replace("B","C"))
        
        # Restrict to just the ranks being considered (ie. when output contains 
        # alternative assignments)
        if "Rank" not in tmp_all.columns:
            tmp_all["Rank"] = 1
            tmp_all["Rel_prob"] = 0
        
        tmp = tmp_all.loc[tmp_all["Rank"].isin(ranks),:].copy()
        tmp["Rank"] = tmp["Rank"].astype(str)
        if "Max_mismatch_m1" in tmp.columns:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob",
                       "Rank","Rel_prob","Dummy_SS","Dummy_res",
                       "Max_mismatch_m1","Max_mismatch_p1",
                       "Num_good_links_m1","Num_good_links_p1"]]
        else:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob",
                       "Rank","Rel_prob","Dummy_SS","Dummy_res"]]
       
        # Convert Res_N column to integer
        tmp.loc[:,"Res_N"] = tmp.loc[:,"Res_N"].fillna(-999)
        tmp.loc[:,"Res_N"] = tmp.loc[:,"Res_N"].astype(int)
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp["Res_type"]=="B"
        tmp.loc[mask,"Res_type"] = "C"
        tmp.loc[mask, "Res_name"] = tmp.loc[mask,"Res_name"].str.replace("B","C")
        # Add a column saying whether a match exists for the spin system 
        # (ie whether a correct assignment is possible)
        tmp["SS_in_pred"] = tmp["SS_name"].isin(tmp_all["Res_name"])
        tmp["Pred_in_SS"] = tmp["Res_name"].isin(tmp_all["SS_name"])
        
        # Add a column with the residue type of the spin system 
        # (as opposed to the prediction assigned to it)
        tmp["SS_type"] = tmp["SS_name"].str[-1:]
        tmp.loc[tmp["Dummy_SS"],"SS_type"] = np.NaN
        
        if assigns is None:
            assigns = tmp
        else:
            assigns = pd.concat([assigns, tmp], ignore_index=True)
        
    # Determine which spin systems were correctly assigned
    assigns["Correct"] = False
    assigns["Status"] = ""
    
    mask = assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), 
                "Correct"] = True
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), 
                "Status"] = "Correctly assigned"
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), 
                "Status"] = "Misassigned"
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Wrongly unassigned"
    
    mask = ~assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), 
                "Status"] = "Wrongly assigned"
    assigns.loc[mask & assigns["Dummy_res"], "Correct"] = True
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Correctly unassigned"
    
    assigns.loc[assigns["Dummy_SS"], "Status"] = "Dummy SS"
    assigns.loc[:,"Status"] = assigns["Status"].astype("category")
    
    # Make the summary dataframe
    ID_unique = assigns["ID"].unique()
    Rank_unique = assigns["Rank"].unique()
    
    # We want:
    # ID1 Rank1
    # ID1 Rank2
    # ...
    # ID1 RankN
    # ID2 Rank1
    # ...
    # ID2 RankN
    # ID3 Rank1
    # ...
    Rank_list = list(Rank_unique)*len(ID_unique)
    ID_list2 = list(ID_unique.repeat(len(Rank_unique)))
    status_list=["Correctly assigned","Correctly unassigned","Dummy SS",
                 "Misassigned","Wrongly assigned","Wrongly unassigned"]
    
    summary = pd.DataFrame({"ID":ID_list2, "Rank":Rank_list}, 
                           columns=["ID", "Rank"]+status_list)

    for i in summary.index:
        tmp = assigns.loc[(assigns["ID"] == summary.loc[i, "ID"]) &
                          (assigns["Rank"] == summary.loc[i, "Rank"]),:]
        for status in status_list:
            summary.loc[i, status] = sum(tmp["Status"]==status)
                
    #summary = summary.loc[:,status_list].fillna(0).astype(int)
    summary["N"] = summary.loc[:,status_list].apply(sum, axis=1)
    summary["N_SS"] = summary["N"] - summary["Dummy SS"]
    
    
    # Add rows with sums for all ranks
    for r in Rank_unique:
        sum_row = (summary.loc[summary["Rank"]==r,:].sum(axis=0).
                       to_frame("Sum").transpose())
        sum_row["ID"] = "Sum"
        sum_row["Rank"] = r
        summary = pd.concat([sum_row, summary], sort=False, ignore_index=True)

    summary["Pc_correct"] = (summary["Correctly assigned"]+
                             summary["Correctly unassigned"]) / summary["N_SS"]
    
    return([assigns, summary])

def collect_assignment_results(data_dir, testset_df, ID_list, prefix="", output_file=None):
    """Function to check assignment accuracy
    
    Returns a list with two items: the complete set of assignments for all 
    proteins, and a summary of the accuracy for each protein
    """
    # Nb. make sure data_dir ends with a forward slash

    assigns = None
    for i in ID_list:
        tmp = pd.read_csv(
                data_dir/(prefix+testset_df.loc[i, "out_name"]+".txt"), 
                sep="\t", index_col=False)
        tmp["ID"] = i
        
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp["Res_type"]=="B"
        tmp.loc[mask,"Res_type"] = "C"
        tmp.loc[mask, "Res_name"] = (tmp.loc[mask,"Res_name"].
                   str.replace("B","C"))
        
        # Tidy up Rank column, and create if missing
        if "Rank" not in tmp.columns:
            tmp["Rank"] = 1
            tmp["Rel_prob"] = 0
        tmp["Rank"] = tmp["Rank"].astype(str)
        
        if "Max_mismatch_m1" in tmp.columns:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob",
                       "Rank","Rel_prob","Dummy_SS","Dummy_res",
                       "Max_mismatch_m1","Max_mismatch_p1",
                       "Num_good_links_m1","Num_good_links_p1", "Confidence"]]
        else:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob",
                       "Rank","Rel_prob","Dummy_SS","Dummy_res"]]
       
        # Convert Res_N column to integer
        tmp.loc[:,"Res_N"] = tmp.loc[:,"Res_N"].fillna(-999)
        tmp.loc[:,"Res_N"] = tmp.loc[:,"Res_N"].astype(int)
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp["Res_type"]=="B"
        tmp.loc[mask,"Res_type"] = "C"
        tmp.loc[mask, "Res_name"] = tmp.loc[mask,"Res_name"].str.replace("B","C")
        # Add a column saying whether a match exists for the spin system 
        # (ie whether a correct assignment is possible)
        tmp["SS_in_pred"] = tmp["SS_name"].isin(tmp["Res_name"])
        tmp["Pred_in_SS"] = tmp["Res_name"].isin(tmp["SS_name"])
        
        # Add a column with the residue type of the spin system 
        # (as opposed to the prediction assigned to it)
        tmp["SS_type"] = tmp["SS_name"].str[-1:]
        tmp.loc[tmp["Dummy_SS"],"SS_type"] = np.NaN
        
        if assigns is None:
            assigns = tmp
        else:
            assigns = pd.concat([assigns, tmp], ignore_index=True)
        
    # Determine which spin systems were correctly assigned
    assigns["Correct"] = False
    assigns["Status"] = ""
    
    mask = assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), 
                "Correct"] = True
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), 
                "Status"] = "Correctly assigned"
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), 
                "Status"] = "Misassigned"
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Wrongly unassigned"
    
    mask = ~assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), 
                "Status"] = "Wrongly assigned"
    assigns.loc[mask & assigns["Dummy_res"], "Correct"] = True
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Correctly unassigned"
    
    assigns.loc[assigns["Dummy_SS"], "Status"] = "Dummy SS"
    assigns.loc[:,"Status"] = assigns["Status"].astype("category")
    
    if output_file is not None:
        assigns.to_csv(output_file, sep="\t", float_format="%.3f")
    
    return(assigns)

def summarise_results(assigns, output_file=None):
    # Make the summary dataframe
    ID_unique = assigns["ID"].unique()
    Rank_unique = assigns["Rank"].unique()
    
    # We want:
    # ID1 Rank1
    # ID1 Rank2
    # ...
    # ID1 RankN
    # ID2 Rank1
    # ...
    # ID2 RankN
    # ID3 Rank1
    # ...
    Rank_list = list(Rank_unique)*len(ID_unique)
    ID_list2 = list(ID_unique.repeat(len(Rank_unique)))
    status_list = ["Correctly assigned","Correctly unassigned","Dummy SS",
                   "Misassigned","Wrongly assigned","Wrongly unassigned"]
    confidence_list = ["High","Medium","Low","Unreliable","Undefined"]
    
    summary = pd.DataFrame({"ID":ID_list2, "Rank":Rank_list}, 
                           columns=["ID", "Rank"]+status_list)

    for i in summary.index:
        tmp = assigns.loc[(assigns["ID"] == summary.loc[i, "ID"]) &
                          (assigns["Rank"] == summary.loc[i, "Rank"]),:]
        for status in status_list:
            summary.loc[i, status] = sum(tmp["Status"]==status)
        for c in confidence_list:
            summary.loc[i, c] = sum(tmp["Confidence"]==c)
            summary.loc[i, c+"_correct"] = sum((tmp["Confidence"]==c) & 
                                               tmp["Correct"])
                
    #summary = summary.loc[:,status_list].fillna(0).astype(int)
    summary["N"] = summary.loc[:,status_list].apply(sum, axis=1)
    summary["N_SS"] = summary["N"] - summary["Dummy SS"]
    
    
    # Add rows with sums for all ranks
    for r in Rank_unique:
        sum_row = (summary.loc[summary["Rank"]==r,:].sum(axis=0).
                       to_frame("Sum").transpose())
        sum_row["ID"] = "Sum"
        sum_row["Rank"] = r
        summary = pd.concat([sum_row, summary], sort=False, ignore_index=True)
    
    summary["Pc_correct"] = 100 * (summary["Correctly assigned"]+
                             summary["Correctly unassigned"]) / summary["N_SS"]
    for c in confidence_list:
        if c=="Dummy_SS":
            continue
        summary[c+"_pc"] = 100 * summary[c]/summary["N_SS"]
        summary[c+"_pc_correct"] = 100 * summary[c+"_correct"]/summary[c].replace({0:np.NaN})
    
    if output_file is not None:
        assigns.to_csv(output_file, sep="\t", float_format="%.3f")
    
    return(summary)

