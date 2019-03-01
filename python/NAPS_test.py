#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Script to test NAPS functionality
Run with a single argument, the path to the NAPS directory
eg "python NAPS_test.py /Users/aph516/GitHub/NAPS/"

@author: aph516
"""

import numpy as np
import pandas as pd
from subprocess import run, Popen, PIPE
from pathlib import Path
import argparse
from plotnine import *

# Set the path to the NAPS directory
parser = argparse.ArgumentParser(
        description="Test script for NAPS (NMR Assignments from Predicted Shifts)")
parser.add_argument("NAPS_path", help="Path to the top-level NAPS directory.")
parser.add_argument("-p", "--python_cmd", default="python")
parser.add_argument("--assign", action="store_true", help="Do NAPS assignment for selected tests.")
parser.add_argument("--analyse", action="store_true", help="Analyse results for selected tests.")
parser.add_argument("-N", default=None, help="Limit to first N datasets.")
parser.add_argument("-t", "--test", nargs="+", default="all", 
                    help="Specify a particular test to run.")

if True:
    args = parser.parse_args()
else:
    #args = parser.parse_args(("C:/Users/kheyam/Documents/GitHub/NAPS/",
    args = parser.parse_args(("/Users/aph516/GitHub/NAPS/",
                              #"--assign",
                              "--analyse",
                              "-t","hadamac",))

path = Path(args.NAPS_path)
#path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")

# Import metadata on the test datasets
testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, 
                           names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/x for x in "data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/x for x in "data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
testset_df.index = testset_df["ID"]

# Define some lists of ID's with particular characteristics
id_all = testset_df["ID"].tolist()  # All ID's

    # Only ID's with info for C, CA and CB (not necessarily complete)
id_all_carbons = [
        'A002', 'A003', 'A004', 'A005', 'A006', 'A008', 'A009', 'A010', 'A011',
        'A012', 'A013', 'A014', 'A015', 'A016', 'A017', 'A018', 'A019', 'A020',
        'A021', 'A023', 'A025', 'A026', 'A027', 'A028', 'A029', 'A033', 'A035',
        'A036', 'A037', 'A039', 'A043', 'A044', 'A045', 'A049', 'A050', 'A051',
        'A053', 'A059', 'A061', 'A062', 'A066', 'A067', 'A069']

id_missing_carbons = list(set(id_all) - set(id_all_carbons))

# Limit how many datasets are assigned/analysed
if args.N is not None:
    id_all = id_all[0:args.N]
    id_all_carbons = id_all_carbons[0:args.N]
#%%

def check_assignment_accuracy(data_dir, ranks=[1], prefix="", ID_list=id_all):
    """Function to check assignment accuracy
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
        if "Max_mismatch_prev" in tmp.columns:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob",
                       "Rank","Rel_prob","Dummy_SS","Dummy_res",
                       "Max_mismatch_prev","Max_mismatch_next",
                       "Num_good_links_prev","Num_good_links_next"]]
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

def save_summary_plot(assigns, summary, out_dir):
    plt = ggplot(data=assigns) #[~assigns["Dummy_SS"]])
    plt = plt + geom_bar(aes(x="ID", fill="Status"), 
                         position=position_fill(reverse=True))
    plt = plt + geom_text(aes(x="summary.index", label="Pc_correct"), y=0.1, 
                          format_string="{:.0%}", data=summary, angle=90)
    plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=0.5))
    #plt = plt + scale_x_discrete(breaks=summary["ID"].tolist())
    plt.save(path/("plots/summary_"+out_dir+".pdf"), 
             height=210, width=297, units="mm")

def save_alt_summary_plots(assigns, summary, out_dir):
    plt = ggplot(data=assigns)
    plt = plt + geom_bar(aes(x="ID", fill="Status"), position=position_fill(reverse=True)) 
    plt = plt + facet_grid("Rank ~ .")
    plt = plt + theme(axis_text_x = element_text(angle=90))
    plt.save(path/"plots"/("summary_"+out_dir+".pdf"), height=210, width=297, units="mm")
    
    tmp = summary
    tmp["Pc_correct"] = tmp["Pc_correct"].astype(float)
    rank_cat = pd.api.types.CategoricalDtype(categories=["3","2","1"], ordered=True)
    tmp["Rank"] = tmp["Rank"].astype(str).astype(rank_cat)
    plt = ggplot(data=tmp) + geom_bar(aes(x="ID", y="Pc_correct", fill="Rank"), stat="identity")
    plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=0.5))
    plt = plt + scale_y_continuous(breaks=np.linspace(0,1,11))
    plt.save(path/"plots"/(out_dir+"_correct.pdf"), height=210, width=297, units="mm")

#%% Test all proteins in the using most basic settings
if "basic" in args.test or "all" in args.test:
    out_dir = "basic"
    if args.assign:
        for i in id_all:
            print((path/("output/testset/"+testset_df.loc[i, "out_name"]+".txt")).as_posix())
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(),
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "-c", (path/"config/config_plot.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix(),
                    "--plot_file", (path/"plots"/out_dir/(testset_df.loc[i, "out_name"]+"_strips.pdf")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_basic, summary_basic = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all)
        summary_basic.to_csv(path/("output/"+out_dir+"_summary.txt") , sep="\t", float_format="%.3f")
    
        save_summary_plot(assigns_basic, summary_basic, out_dir)
        
        summary_basic.iloc[1:,:]["Pc_correct"].quantile([0,0.25,0.5,0.75,1])
        
        # Make a matrix showing how often each residue type is misassigned to a different type
        misassigned_types_basic = assigns_basic[["SS_type","Res_type"]].groupby(["SS_type","Res_type"]).size().unstack().fillna(0)
        type_count = misassigned_types_basic.sum(axis=1)    # Count occurrences of each residue type in observations
        np.fill_diagonal(misassigned_types_basic.values, 0) # Set diagonal to zero
        # Calculate % of each type that gets misassigned to each other type
        misassigned_types2_basic = (misassigned_types_basic/type_count)*100
        
        tmp = assigns_basic[assigns_basic["Dummy_SS"]==False].groupby("Res_type")["Correct"]
        tmp.sum()/tmp.count()    
        
        tmp = assigns_basic[(assigns_basic["Dummy_SS"]==False) & (assigns_basic["Status"]=="Misassigned")]
        tmp["Type_match"] = tmp["SS_type"]==tmp["Res_type"]
        tmp2 = tmp.groupby("SS_type")["Type_match"]
        tmp_basic=(tmp2.sum()/tmp2.count()).sort_values(ascending=False)
        
        tmp3 = tmp[(tmp["Type_match"]==False) & (tmp["Status"]=="Misassigned")]
        tmp3.groupby("SS_type")["ID"].count()
        
        # Check if there's any pattern in the log_probabilities
        tmp = assigns_basic[(assigns_basic["ID"]=="A006") & 
                          ~assigns_basic["Dummy_SS"] & 
                          ~assigns_basic["Dummy_res"]]
        (ggplot(tmp) + geom_density(aes(x="Log_prob", colour="Correct")) 
        + xlim(-100, 0) )
        
        tmp = assigns_basic[~assigns_basic["Dummy_SS"] & 
                            ~assigns_basic["Dummy_res"]]
        (ggplot(tmp) + geom_boxplot(aes(y="Log_prob", x="ID", colour="Correct")) 
        + ylim(-100,0) )
    
    
#%% Test effect of correcting the predicted shifts
if "pred_correction" in args.test or "all" in args.test:
    out_dir = "pred_correction"
    if args.assign:
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(),
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "-c", (path/"config/config_pred_correction.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_pc, summary_pc = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all)
        summary_pc.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_pc, summary_pc, out_dir)


#%% Test effect of accounting for correlated errors
if "delta_correlation" in args.test or "all" in args.test:
    out_dir = "delta_correlation"
    if args.assign:
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(),
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "-c", (path/"config/config_delta_corr.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_dc, summary_dc = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all)
        summary_dc.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_dc, summary_dc, out_dir)

#%% Test effect of accounting for correlated errors *and* correcting the predicted shifts
if "delta_correlation2" in args.test or "all" in args.test:
    out_dir = "delta_correlation2"
    if args.assign:
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(),
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "-c", (path/"config/config_delta_corr2.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_dc2, summary_dc2 = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all)
        summary_dc2.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_dc2, summary_dc2, out_dir)
        

#%% Test effect of including HADAMAC amino acid type information
if "hadamac" in args.test or "all" in args.test:
    out_dir = "hadamac"
    if args.assign:
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(),
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML",
                    "-c", (path/"config/config_hadamac.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_hadamac, summary_hadamac = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all)
        summary_hadamac.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_hadamac, summary_hadamac, out_dir)
       
        summary_hadamac.iloc[1:,:]["Pc_correct"].quantile([0,0.25,0.5,0.75,1])
        
        tmp = assigns_hadamac[(assigns_hadamac["Dummy_SS"]==False) & (assigns_hadamac["Status"]=="Misassigned")]
        tmp["Type_match"] = tmp["SS_type"]==tmp["Res_type"]
        tmp2 = tmp.groupby("SS_type")["Type_match"]
        tmp_had=(tmp2.sum()/tmp2.count()).sort_values(ascending=False)
        
        # Make a matrix showing how often each residue type is misassigned to a different type
        misassigned_types_hadamac = assigns_hadamac[["SS_type","Res_type"]].groupby(["SS_type","Res_type"]).size().unstack().fillna(0)
        type_count_hadamac = misassigned_types_hadamac.sum(axis=1)    # Count occurrences of each residue type in observations
        np.fill_diagonal(misassigned_types_hadamac.values, 0) # Set diagonal to zero
        # Calculate % of each type that gets misassigned to each other type
        misassigned_types2_hadamac = (misassigned_types_hadamac/type_count_hadamac)*100
        
        
#%% Test alternative assignments with and without HADAMAC
if "alt_assign" in args.test or "all" in args.test:
    out_dir = "alt_assign"
    if args.assign:
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "-c", (path/"config/config_alt_assign.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_alt, summary_alt = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt = assigns_alt.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)

if "alt_hadamac" in args.test or "all" in args.test:
    out_dir = "alt_hadamac"
    if args.assign:
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2",
                    "--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML",
                    "-c", (path/"config/config_alt_hadamac.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:        
        assigns_alt, summary_alt = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt = assigns_alt.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_hadamac, summary_hadamac, out_dir)
    
#%% Test alternative assignments with reduced atom types
#### HNCO with and without HADAMAC
if "alt_hnco" in args.test or "all" in args.test:
    out_dir = "alt_hnco"
    if args.assign:
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2", 
                    "-c", (path/"config/config_alt_hnco.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
            
    if args.analyse:
        assigns_alt_hnco, summary_alt_hnco = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt_hnco = assigns_alt_hnco.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt_hnco, summary_alt_hnco, out_dir)

if "alt_hnco2" in args.test or "all" in args.test:
    out_dir = "alt_hnco2"
    if args.assign:
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2", 
                    "-c", (path/"config/config_alt_hnco2.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
            
    if args.analyse:
        assigns_alt_hnco2, summary_alt_hnco2 = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt_hnco2 = assigns_alt_hnco.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco2.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco2.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt_hnco2, summary_alt_hnco2, out_dir)

#### HNCO + HNCACB
if "alt_hnco_hncacb" in args.test or "all" in args.test:
    out_dir = "alt_hnco_hncacb"
    if args.assign:        
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2", 
                    "-c", (path/"config/config_alt_hnco_hncacb.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
            
    if args.analyse:        
        assigns_alt_hnco_hncacb, summary_alt_hnco_hncacb = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt_hnco_hncacb = assigns_alt_hnco_hncacb.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco_hncacb.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco_hncacb.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)

# All CA and CO shifts, but no CB
if "alt_ca_co" in args.test or "all" in args.test:
    out_dir = "alt_ca_co"
    if args.assign:    
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
                    testset_df.loc[i, "obs_file"].as_posix(), 
                    testset_df.loc[i, "preds_file"].as_posix(), 
                    (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                    "--shift_type", "test",
                    "--pred_type", "shiftx2", 
                    "-c", (path/"config/config_alt_ca_co.txt").as_posix(),
                    "-l", (path/"output"/out_dir/(testset_df.loc[i, "out_name"]+".log")).as_posix()]
            run(cmd)
    
    if args.analyse:            
        assigns_alt_ca_co, summary_alt_ca_co = check_assignment_accuracy(path/"output"/out_dir, ID_list=id_all_carbons)
        assigns_alt_ca_co = assigns_alt_ca_co.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_ca_co.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_ca_co.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)
