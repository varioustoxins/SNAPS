#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Script to test NAPS functionality
Run with a single argument, the path to the NAPS directory
eg "python NAPS_test.py /Users/aph516/Github/NAPS/"

@author: aph516
"""

import pandas as pd
from subprocess import run, Popen, PIPE
from pathlib import Path
import argparse
from plotnine import *

# Set the path to the NAPS directory
parser = argparse.ArgumentParser(description="Test script for NAPS (NMR Assignments from Predicted Shifts)")
parser.add_argument("NAPS_path", help="Path to the top-level NAPS directory.")
parser.add_argument("-p", "--python_cmd", default="python")
parser.add_argument("-t", "--test", default="all", help="Specify a particular test to run.")
parser.add_argument("-N", "--N_tests", default=61, type=int, help="only process the first N datasets")

if False:
    #args = parser.parse_args("C:/Users/kheyam/Documents/GitHub/NAPS/ -N 10".split())
    args = parser.parse_args("/Users/aph516/GitHub/NAPS/ -t alt_assignments".split())
else:
    args = parser.parse_args()

path = Path(args.NAPS_path)
#path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")

# Import metadata on the test datasets
testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/x for x in "data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/x for x in "data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
testset_df.index = testset_df["ID"]

#%%

def check_assignment_accuracy(data_dir, ranks=[1], prefix="", N=61):
    """Function to check assignment accuracy
    """
    # Nb. make sure data_dir ends with a forward slash

    assigns = None
    for i in testset_df.index[0:N]:
        tmp_all = pd.read_csv(data_dir/(prefix+testset_df.loc[i, "out_name"]+".txt"), sep="\t", index_col=0)
        tmp_all["ID"] = i
        
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp_all["Res_type"]=="B"
        tmp_all.loc[mask,"Res_type"] = "C"
        tmp_all.loc[mask, "Res_name"] = tmp_all.loc[mask,"Res_name"].str.replace("B","C")
        
        # Restrict to just the ranks being considered (ie. when output contains alternative assignments)
        if "Rank" not in tmp_all.columns:
            tmp_all["Rank"] = 1
            tmp_all["Rel_prob"] = 0
        
        tmp = tmp_all.loc[tmp_all["Rank"].isin(ranks),:].copy()
        tmp["Rank"] = tmp["Rank"].astype(str)
        if "Max_mismatch_prev" in tmp.columns:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob","Rank","Rel_prob","Dummy_SS","Dummy_res","Max_mismatch_prev","Max_mismatch_next","Num_good_links_prev","Num_good_links_next"]]
        else:
            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob","Rank","Rel_prob","Dummy_SS","Dummy_res"]]
#        else:
#            tmp = tmp_all
#            tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Log_prob","Dummy_SS","Dummy_res"]]
        
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
        
        if type(assigns)!=type(tmp):
            assigns = tmp
        else:
            assigns = pd.concat([assigns, tmp], ignore_index=True)
    
    # Determine which spin systems were correctly assigned
    assigns["Correct"] = False
    assigns["Status"] = ""
    
    mask = assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), "Correct"] = True
    assigns.loc[mask & (assigns["SS_name"]==assigns["Res_name"]), "Status"] = "Correctly assigned"
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), "Status"] = "Misassigned"
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Wrongly unassigned"
    
    mask = ~assigns["SS_in_pred"] & ~assigns["Dummy_SS"]
    assigns.loc[mask & (assigns["SS_name"]!=assigns["Res_name"]), "Status"] = "Wrongly assigned"
    assigns.loc[mask & assigns["Dummy_res"], "Correct"] = True
    assigns.loc[mask & assigns["Dummy_res"], "Status"] = "Correctly unassigned"
    
    assigns.loc[assigns["Dummy_SS"], "Status"] = "Dummy SS"
    assigns.loc[:,"Status"] = assigns["Status"].astype("category")
    
    #assigns.groupby(["ID","Status"])["Status"].count()
#    if "Rank" in tmp_all.columns:
    
    summary = pd.DataFrame(columns=["Rank","Correctly assigned","Correctly unassigned","Dummy SS","Misassigned","Wrongly assigned","Wrongly unassigned","N","N_SS","Pc_correct"])
    summary = summary.append(assigns.groupby(["ID","Rank","Status"])["Status"].count().unstack(2), sort=False)
#    else:
#    summary = pd.DataFrame(columns=["Correctly assigned","Correctly unassigned","Dummy SS","Misassigned","Wrongly assigned","Wrongly unassigned","N","N_SS","Pc_correct"])
#    summary = summary.append(assigns.groupby(["ID","Status"])["Status"].count().unstack(1), sort=False)
        
    summary = summary.fillna(0).astype(int)
    summary["N"] = summary.apply(sum, axis=1)
    summary["N_SS"] = summary["N"] - summary["Dummy SS"]
    summary = pd.concat([summary.sum(axis=0).to_frame("Sum").transpose(), summary], sort=False)
    summary["Pc_correct"] = (summary["Correctly assigned"]+summary["Correctly unassigned"]) / summary["N_SS"]
    
    return([assigns, summary])

#%%
#### Test all proteins in the using most basic settings
if args.test in ("basic", "all"):
    for i in testset_df.index[0:args.N_tests]:
        print((path/("output/testset/"+testset_df.loc[i, "out_name"]+".txt")).as_posix())
        cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(), 
                testset_df.loc[i, "obs_file"].as_posix(), 
                testset_df.loc[i, "preds_file"].as_posix(), 
                (path/("output/testset/"+testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                "-c", (path/"config/config_plot.txt").as_posix(),
                "-l", (path/("output/testset/"+testset_df.loc[i, "out_name"]+".log")).as_posix(),
                "-p", (path/("plots/testset/"+testset_df.loc[i, "out_name"]+"_strips.pdf")).as_posix()]
        run(cmd)
        
assigns_std, summary_std = check_assignment_accuracy(path/"output/testset/", N=args.N_tests)
summary_std.to_csv(path/"output/testset_summary.txt", sep="\t", float_format="%.3f")

plt = ggplot(data=assigns_std) + geom_bar(aes(x="ID", fill="Status"), position=position_fill(reverse=True))
plt = plt + geom_text(aes(x="summary_std.index", label="Pc_correct"), y=0.1, format_string="{:.0%}", data=summary_std, angle=90)
plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=1))
plt.save(path/"plots/testset_summary.pdf", height=210, width=297, units="mm")
#%%
#### Test effect of accounting for correlated errors
if args.test in ("delta_correlation", "all"):
    for i in testset_df.index[0:args.N_tests]:
        print(testset_df.loc[i, "out_name"])    
        cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(), 
                testset_df.loc[i, "obs_file"].as_posix(), 
                testset_df.loc[i, "preds_file"].as_posix(), 
                (path/("output/delta_correlation/"+testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                "-c", (path/"config/config_delta_corr.txt").as_posix(),
                "-l", (path/("output/delta_correlation/"+testset_df.loc[i, "out_name"]+".log")).as_posix()]
        run(cmd)
        
    assigns_dc, summary_dc = check_assignment_accuracy(path/"output/delta_correlation/", N=args.N_tests)
    summary_dc.to_csv(path/"output/delta_correlation_summary.txt", sep="\t", float_format="%.3f")
    
    plt = ggplot(data=assigns_dc) + geom_bar(aes(x="ID", fill="Status"), position=position_fill(reverse=True))
    plt = plt + geom_text(aes(x="summary_dc.index", label="Pc_correct"), y=0.1, format_string="{:.0%}", data=summary_std, angle=90)
    plt = plt + theme(axis_text_x=element_text(rotation=90, hjust=1))
    plt.save(path/"plots/delta_correlation_summary.pdf", height=210, width=297, units="mm")
#%%
#### Test alternative assignments
if args.test in ("alt_assignments", "all"):
    for i in testset_df.index[0:args.N_tests]:
        print(testset_df.loc[i, "out_name"])    
        cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(), 
                testset_df.loc[i, "obs_file"].as_posix(), 
                testset_df.loc[i, "preds_file"].as_posix(), 
                (path/("output/alt_assign/"+testset_df.loc[i, "out_name"]+".txt")).as_posix(), 
                "-c", (path/"config/config_alt_assign.txt").as_posix(),
                "-l", (path/("output/alt_assign/"+testset_df.loc[i, "out_name"]+".log")).as_posix()]
        run(cmd)
        
    assigns_alt, summary_alt = check_assignment_accuracy(path/"output/alt_assign/", ranks=[1,2,3], N=args.N_tests)
#    assigns_alt["Rank"] = 1
#    summary_alt["Rank"] = 1
#    for i in (2,3):
#        print(i)
#        tmp1, tmp2 = check_assignment_accuracy(path/"output/alt_assign/", ranks=[i], N=args.N_tests)
#        tmp1["Rank"] = i
#        tmp2["Rank"] = i
#        print(tmp1.iloc[0:10, 0:5])
#        assigns_alt = pd.concat([assigns_alt, tmp1], ignore_index=True)
#        summary_alt = pd.concat([summary_alt,tmp2], ignore_index=True)
    assigns_alt = assigns_alt.sort_values(by=["ID", "SS_name", "Rank"])
    
    assigns_alt.to_csv(path/"output/alt_assign_all.txt", sep="\t", float_format="%.3f")
    summary_alt.to_csv(path/"output/alt_assign_summary.txt", sep="\t", float_format="%.3f")
    
    plt = ggplot(data=assigns_alt) + geom_bar(aes(x="ID", fill="Status"), position=position_fill(reverse=True)) + facet_grid("Rank ~ .")
    plt.save(path/"plots/alt_assign_summary.pdf", height=210, width=297, units="mm")
    
    
    tmp1, tmp2 = check_assignment_accuracy(path/"output/alt_assign/", ranks=[1,2,3])
    