#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Script to test NAPS functionality
Run with a single argument, the path to the NAPS directory
eg "python NAPS_test.py /Users/aph516/GitHub/NAPS/"

/anaconda3/bin/python3 NAPS_test.py .. -p /anaconda3/bin/python3 --assign -t basic

@author: aph516
"""

import numpy as np
import pandas as pd
from NAPS_analyse import import_testset_metadata, check_assignment_accuracy, collect_assignment_results, summarise_results
from subprocess import run
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
    args = parser.parse_args(("..",
                              #"--assign",
                              "--analyse"))

path = Path(args.NAPS_path)
#path = Path("..")
#path = Path("/Users/aph516/GitHub/NAPS/")
#path = Path("C:/Users/kheyam/Documents/GitHub/NAPS/")

# Import metadata on the test datasets
testset_df = import_testset_metadata(path)

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
    id_all = id_all[0:int(args.N)]
    id_all_carbons = id_all_carbons[0:int(args.N)]
#%% Some functions

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
    return(plt)

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
    return(plt)
    
def make_cmd(id, out_dir, config_file="config.txt", extra_args=[]):
    """Extra_args is a list of additional command line arguments"""
    cmd = [args.python_cmd, (path/"python/NAPS.py").as_posix(),
            testset_df.loc[id, "obs_file"].as_posix(), 
            testset_df.loc[id, "preds_file"].as_posix(),
            (path/"output"/out_dir/(testset_df.loc[id, "out_name"]+".txt")).as_posix(),
            "--shift_type", "test",
            "--pred_type", "shiftx2",
            "-c", (path/"config"/config_file).as_posix(),
            "-l", (path/"output"/out_dir/(testset_df.loc[id, "out_name"]+".log")).as_posix()]
    cmd = cmd + extra_args
            
    return(cmd)


#%% Test all proteins in the using most basic settings
if "basic" in args.test or "all" in args.test:
    out_dir = "basic"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print((path/("output/testset/"+testset_df.loc[i, "out_name"]+".txt")).as_posix())
            cmd = make_cmd(i, out_dir, "config_plot.txt",
                           ["--plot_file", (path/"plots"/out_dir/(testset_df.loc[i, "out_name"]+"_strips.html")).as_posix()])
            run(cmd)
    
    if args.analyse:        
        #assigns_basic, summary_basic = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
        assigns_basic = collect_assignment_results(path/"output"/out_dir, testset_df, ID_list=id_all)
        summary_basic = summarise_results(assigns_basic)
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
        
        #### Make figures for poster
        # Figure breaking down assignment accuracy
        tmp = assigns_basic[~assigns_basic["Dummy_SS"]]
        plt = ggplot(tmp) + geom_bar(aes(x="Status", fill="Status",
                                     y="100*stat(count)/"+str(len(tmp.index))))
        plt += xlab("Assignment status") + ylab("Percentage")
        plt += scale_fill_brewer("qual", palette=6)
        plt += scale_y_continuous(breaks=np.arange(0,100,10))
        plt += theme_bw() + theme(axis_text_x = element_text(angle=90))
        plt.save(path/"plots/Poster basic accuracy.pdf", height=100, width=100, units="mm")
        
        # Figure showing accuracy vs number of good seq links
        # Note that NA values are always false in > or < comparisons.        
        plt = ggplot(tmp) + geom_bar(aes(x="((Num_good_links_prev>=2) & (Max_mismatch_prev<0.1)).astype(int)"+
                                         "+((Num_good_links_next>=2) & (Max_mismatch_next<0.1)).astype(int)", 
                                     y="100*stat(count)/"+str(len(tmp.index)), 
                                     fill="Correct"))
        plt += xlab("Number of good seq links") + ylab("Percentage")
        plt += scale_fill_brewer("qual", palette=6)
        plt += theme_bw()
        plt.save(path/"plots/Poster basic seq links.pdf", height=100, width=100, units="mm")
        
        # Some alternate ways of assessing how good the seq links are
#        "(Num_good_links_prev==3).astype(int)+(Num_good_links_next==3).astype(int)"
#        ("(Max_mismatch_prev<0.1).astype(int)"+
#         "+(Max_mismatch_next<0.1).astype(int)")
#        ("((Num_good_links_prev>=2) & (Max_mismatch_prev<0.1)).astype(int)"+
#         "+((Num_good_links_next>=2) & (Max_mismatch_next<0.1)).astype(int)")
        
        # Figure showing distribution of log_probabilities
        plt = ggplot(tmp[~tmp["Dummy_res"]])
        plt += geom_density(aes(x="Log_prob", fill="Correct"), alpha=0.5)
        plt += xlim(-50,0)
        plt += scale_fill_brewer(type="qual", palette=6)
        plt += theme_bw()
        plt.save(path/"plots/Poster log_prob distribution.pdf", height=100, width=100, units="mm")
        
#%% Test effect of correcting the predicted shifts
if "pred_correction" in args.test or "all" in args.test:
    out_dir = "pred_correction"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_pred_correction.txt")
            run(cmd)
    
    if args.analyse:        
        assigns_pc, summary_pc = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
        summary_pc.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_pc, summary_pc, out_dir)


#%% Test effect of accounting for correlated errors
if "delta_correlation" in args.test or "all" in args.test:
    out_dir = "delta_correlation"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_delta_corr.txt")
            run(cmd)
    
    if args.analyse:        
        assigns_dc, summary_dc = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
        summary_dc.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_dc, summary_dc, out_dir)

#%% Test effect of accounting for correlated errors *and* correcting the predicted shifts
if "delta_correlation2" in args.test or "all" in args.test:
    out_dir = "delta_correlation2"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_delta_corr2.txt")
            run(cmd)
    
    if args.analyse:        
        assigns_dc2, summary_dc2 = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
        summary_dc2.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_dc2, summary_dc2, out_dir)
        

#%% Test effect of including HADAMAC amino acid type information
if "hadamac" in args.test or "all" in args.test:
    out_dir = "hadamac"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_hadamac.txt",
                           ["--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML"])
            run(cmd)
    
    if args.analyse:        
        assigns_hadamac, summary_hadamac = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
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
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_assign.txt")
            run(cmd)
    
    if args.analyse:        
        assigns_alt, summary_alt = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt = assigns_alt.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)

if "alt_hadamac" in args.test or "all" in args.test:
    out_dir = "alt_hadamac"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_hadamac.txt",
                           ["--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML"])
            run(cmd)
    
    if args.analyse:        
        assigns_alt_hadamac, summary_alt_hadamac = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt_hadamac = assigns_alt_hadamac.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hadamac.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hadamac.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_hadamac, summary_hadamac, out_dir)
    
#%% Test alternative assignments with reduced atom types
#### HNCO with and without HADAMAC
if "alt_hnco" in args.test or "all" in args.test:
    out_dir = "alt_hnco"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_hnco.txt",
                           ["--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML"])
            run(cmd)
            
    if args.analyse:
        assigns_alt_hnco, summary_alt_hnco = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt_hnco = assigns_alt_hnco.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt_hnco, summary_alt_hnco, out_dir)

if "alt_hnco2" in args.test or "all" in args.test:
    out_dir = "alt_hnco2"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_hnco2.txt",
                           ["--test_aa_classes", "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML"])
            run(cmd)
            
    if args.analyse:
        assigns_alt_hnco2, summary_alt_hnco2 = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt_hnco2 = assigns_alt_hnco2.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco2.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco2.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt_hnco2, summary_alt_hnco2, out_dir)

#### HNCO + HNCACB
if "alt_hnco_hncacb" in args.test or "all" in args.test:
    out_dir = "alt_hnco_hncacb"
    if args.assign:    
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_hnco_hncacb.txt")
            run(cmd)
            
    if args.analyse:        
        assigns_alt_hnco_hncacb, summary_alt_hnco_hncacb = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt_hnco_hncacb = assigns_alt_hnco_hncacb.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_hnco_hncacb.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_hnco_hncacb.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)

# All CA and CO shifts, but no CB
if "alt_ca_co" in args.test or "all" in args.test:
    out_dir = "alt_ca_co"
    if args.assign:    
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all_carbons:
            print(testset_df.loc[i, "out_name"])    
            cmd = make_cmd(i, out_dir, "config_alt_ca_co.txt")
            run(cmd)
    
    if args.analyse:            
        assigns_alt_ca_co, summary_alt_ca_co = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all_carbons, ranks=[1,2,3])
        assigns_alt_ca_co = assigns_alt_ca_co.sort_values(by=["ID", "SS_name", "Rank"])
        
        assigns_alt_ca_co.to_csv(path/("output/"+out_dir+"_all.txt"), sep="\t", float_format="%.3f")
        summary_alt_ca_co.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
    
        save_alt_summary_plots(assigns_alt, summary_alt, out_dir)

#%% Test iterated assignment
if "iterated" in args.test or "all" in args.test:
    out_dir = "iterated"
    if args.assign:
        # Create output directory, if it doesn't already exist
        (path/"output"/out_dir).mkdir(parents=True, exist_ok=True)
        for i in id_all:
            print((path/("output/testset/"+testset_df.loc[i, "out_name"]+".txt")).as_posix())
            cmd = make_cmd(i, out_dir, "config.txt", ["--iterated"])
            run(cmd)
    
    if args.analyse:        
        assigns_iter, summary_iter = check_assignment_accuracy(path/"output"/out_dir, testset_df, ID_list=id_all)
        summary_iter.to_csv(path/("output/"+out_dir+"_summary.txt"), sep="\t", float_format="%.3f")
        
        save_summary_plot(assigns_iter, summary_iter, out_dir)

#%% Test stuff

# Try to work out why A063 and some others are so inaccurate
#tmp = summary_basic[["ID","N","Dummy SS"]]
#tmp["Dummy_res"] = summary_basic["Correctly unassigned"] + summary_basic["Wrongly unassigned"]
#tmp["basic"] = summary_basic["Pc_correct"]
#tmp["iter"] = summary_iter["Pc_correct"]
#tmp["diff"] = tmp["iter"] - tmp["basic"]
#
#obs_only = {}
#preds_only = {}
#for i in assigns_basic["ID"].unique():
#    Res_list = assigns_basic.loc[assigns_basic["ID"]==i, "Res_name"]
#    SS_list = assigns_basic.loc[assigns_basic["ID"]==i, "SS_name"]
#    preds_only[i] = len(set(Res_list).difference(SS_list))
#    obs_only[i] = len(set(SS_list).difference(Res_list))
#
#diffs = pd.DataFrame({"ID":list(obs_only.keys()), "obs_only":list(obs_only.values()),
#                      "preds_only":list(preds_only.values())})
#
#tmp = tmp.merge(diffs, how="outer")
#
#ggplot(tmp) + geom_point(aes(x="N_diffs", y="iter.astype(float)"))
