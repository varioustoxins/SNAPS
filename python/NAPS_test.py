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

# Set the path to the NAPS directory
parser = argparse.ArgumentParser(description="Test script for NAPS (NMR Assignments from Predicted Shifts)")
parser.add_argument("NAPS_path", help="Path to the top-level NAPS directory.")
parser.add_argument("-p", "--python_cmd", default="python")
parser.add_argument("-t", "--test", default="all", help="Specify a particular test to run.")
parser.add_argument("-N", "--N_tests", default=61, type=int, help="only process the first N datasets")
args = parser.parse_args()
if False:
    args = parser.parse_args("C:/Users/kheyam/Documents/GitHub/NAPS/ -N 10".split())
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
        tmp = pd.read_csv(data_dir/(prefix+testset_df.loc[i, "out_name"]+".txt"), sep="\t", index_col=0)
        if "Rank" in tmp.columns:
            tmp = tmp.loc[tmp["Rank"].isin(ranks),:]
        tmp["ID"] = i
        tmp = tmp[["ID","Res_N","Res_type","Res_name","SS_name","Dummy_SS","Dummy_res"]]
        # Convert Res_N column to integer
        tmp["Res_N"] = tmp["Res_N"].fillna(-999)
        tmp["Res_N"] = tmp["Res_N"].astype(int)
        # Cysteines in disulphide bridges are often named B in this dataset. 
        # We don't need to know this, so change to C
        mask = tmp["Res_type"]=="B"
        tmp.loc[mask,"Res_type"] = "C"
        tmp.loc[mask, "Res_name"] = tmp.loc[mask,"Res_N"].astype(str)+tmp.loc[mask, "Res_type"]
        # Add a column saying whether a match exists for the spin system 
        # (ie whether a correct assignment is possible)
        tmp["SS_in_pred"] = tmp["SS_name"].isin(tmp["Res_name"])
        tmp["Pred_in_SS"] = tmp["Res_name"].isin(tmp["SS_name"])
        
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
    
    #assigns.groupby(["ID","Status"])["Status"].count()
    summary = pd.DataFrame(columns=["Correctly assigned","Correctly unassigned","Dummy SS","Misassigned","Wrongly assigned","Wrongly unassigned","N","N_SS","Pc_correct"])
    summary = summary.append(assigns.groupby(["ID","Status"])["Status"].count().unstack(1), sort=False)
        
    summary = summary.fillna(0).astype(int)
    summary["N"] = summary.apply(sum, axis=1)
    summary["N_SS"] = summary["N"] - summary["Dummy SS"]
    summary = pd.concat([summary.sum(axis=0).to_frame("Sum").transpose(), summary])
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
                "-c", (path/"config/config.txt").as_posix(),
                "-l", (path/("output/testset/"+testset_df.loc[i, "out_name"]+".log")).as_posix()]
        run(cmd)
        
assigns, summary = check_assignment_accuracy(path/"output/testset/", N=args.N_tests)
summary.to_csv(path/"output/testset_summary.txt", sep="\t", float_format="%.3f")
#%%
#### Test effect of accounting for correlated errors
if args.test in ("delta_correlation", "all"):
    for i in testset_df.index[0:args.N_tests]:
        print(testset_df.loc[i, "out_name"])    
        cmd = [args.python_cmd, args.python_cmd, (path/"python/NAPS.py").as_posix(), 
                testset_df.loc[i, "obs_file"].as_posix(), 
                testset_df.loc[i, "preds_file"].as_posix(), 
                (path/("output/delta_correlation/"+testset_df.loc[i, "out_name"]+".txt")).as_posix(),
                "-c", (path/"config/config.txt").as_posix(),
                "-l", (path/("output/delta_correlation/"+testset_df.loc[i, "out_name"]+".log")).as_posix(),
                "--delta_correlation"]
        run(cmd)
        
    assigns2, summary2 = check_assignment_accuracy(path/"output/delta_correlation/", N=args.N_tests)
    summary2.to_csv(path/"output/delta_correlation_summary.txt", sep="\t", float_format="%.3f")

#%%
#### Test alternative assignments
if args.test in ("alt_assignments", "all"):
    for i in testset_df.index[0:args.N_tests]:
        print(testset_df.loc[i, "out_name"])    
        cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(), 
                testset_df.loc[i, "obs_file"].as_posix(), 
                testset_df.loc[i, "preds_file"].as_posix(), 
                (path/("output/testset/alt_"+testset_df.loc[i, "out_name"]+".txt")).as_posix(), 
                "-c", (path/"config/config.txt").as_posix(),
                "-l", (path/("output/testset/alt_"+testset_df.loc[i, "out_name"]+".log")).as_posix(),
                "--delta_correlation", "--alt_assignments", "2"]
        run(cmd)
        
    assigns3, summary3 = check_assignment_accuracy(path/"output/testset/", prefix="alt_", ranks=[2], N=args.N_tests)
    summary3.to_csv(path/"output/alt_assign_summary.txt", sep="\t", float_format="%.3f")