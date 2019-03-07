# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 21:32:56 2019

@author: Alex
"""

import pandas as pd
import numpy as np
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
from scipy.optimize import linear_sum_assignment
import argparse
from pathlib import Path
from distutils.util import strtobool
import logging

#%% Command line arguments

parser = argparse.ArgumentParser(description="NAPS (NMR Assignments from Predicted Shifts)")
parser.add_argument("input_file", 
                    help="""A table of observed chemical shifts.""")
parser.add_argument("output_file")
parser.add_argument("--shift_type", 
                    choices=["naps", "ccpn", "sparky", 
                             "xeasy", "nmrpipe", "test"], 
                    default="naps", 
                    help="""The format of the observed shift file.""")
parser.add_argument("--pred_file")
parser.add_argument("--pred_type", 
                    choices=["shiftx2", "sparta+"],
                    default="shiftx2", 
                    help="The file format of the predicted shifts")
parser.add_argument("-c", "--config_file", 
                    default="/Users/aph516/GitHub/NAPS/python/config.txt")
parser.add_argument("-l", "--log_file", default=None)
parser.add_argument("--plot_file") 


#args = parser.parse_args()
args = parser.parse_args(("../data/testset/simplified_BMRB/4834.txt "+
                "../output/NAPS_iter_test.txt " +
                "--shift_type test " +
                "--pred_file ../data/testset/shiftx2_results/A003_1LM4B.cs "+
                "--pred_type shiftx2 "+
                "-c ../config/config.txt "+
                "-l ../output/NAPS_iter_test.log "+
                "--plot_file ../plots/NAPS_iter_test.pdf").split())

#%% Set up logging
if isinstance(args.log_file, str):
    print(args.log_file)
    logging.basicConfig(filename=args.log_file, filemode='w', level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s")
else:
    logging.basicConfig(level=logging.ERROR)

#%% Main program
    
a = NAPS_assigner()

# Import config file
a.read_config_file(args.config_file)
logging.info("Read in configuration from %s.", args.config_file)

# Import observed shifts
importer = NAPS_importer()

if args.shift_type=="test":
    importer.import_testset_shifts(args.input_file)
else:
    importer.import_obs_shifts(args.input_file, args.shift_type, SS_num=False)
a.obs = importer.obs
logging.info("Read in %d spin systems from %s.", 
             len(a.obs["SS_name"]), args.input_file)

# Import predicted shifts
a.import_pred_shifts(args.pred_file, args.pred_type)
logging.info("Read in %d predicted residues from %s.", 
             len(a.preds["Res_name"]), args.pred_file)

# Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix2(sf=1, verbose=False)
logging.info("Calculated log probability matrix (%dx%d).", 
             a.log_prob_matrix.shape[0], a.log_prob_matrix.shape[1])
matching1 = a.find_best_assignments()
logging.info("Calculated best assignment.")
assign_df = a.make_assign_df(matching1, set_assign_df=True)
assign_df = a.check_assignment_consistency(threshold=0.1)
logging.info("Checked assignment consistency.")

obs = a.obs
preds = a.preds

#%% 1st attempt: try excluding each inconsistent match in turn, and keep going if consistency improves.

last_consistent=-1
max_consistent = ((assign_df["Num_good_links_prev"]==3) & 
                  (assign_df["Num_good_links_next"]==3)).sum()
count=1

while max_consistent>last_consistent:
    print("Round %d" % count)
    count = count+1
    last_consistent = max_consistent
    # Work out which residues are consistent
    consistent = a.assign_df.loc[(a.assign_df["Num_good_links_prev"]==3) & 
                               (a.assign_df["Num_good_links_next"]==3),
                               ["SS_name", "Res_name"]]
    inconsistent = a.assign_df.loc[~a.assign_df["Res_name"].isin(consistent["Res_name"]), ["SS_name", "Res_name"]]
    
    alt_assigns = pd.DataFrame({"SS_name":inconsistent["SS_name"],
                                "Res_name":inconsistent["Res_name"],
                                "Sum_prob":np.NaN, "Total_consistent":np.NaN}, 
                                index=inconsistent.index)
    matching_dict={}
    
    
    for i in inconsistent.index:
        matching2 = a.find_best_assignment2(inc=consistent, exc=inconsistent.loc[[i],:])
        matching_dict[i] = matching2
        alt_assigns.loc[i,"Sum_prob"] = sum(a.log_prob_matrix.lookup(matching2["SS_name"], matching2["Res_name"]))
        tmp = a.make_assign_df(matching2, set_assign_df=False)
        tmp = a.check_assignment_consistency(assign_df=tmp)
        alt_assigns.loc[i,"Total_consistent"] = ((tmp["Num_good_links_prev"]==3) & 
                                                 (tmp["Num_good_links_next"]==3)).sum()
    
    max_consistent = alt_assigns["Total_consistent"].max()
    mask = alt_assigns["Total_consistent"]==max_consistent
    i = alt_assigns.loc[mask,"Sum_prob"].idxmax()
    print("Top result: %s\tTotal consistent: %d, Sum_probability: %d" % 
          (alt_assigns.loc[i, "Res_name"],alt_assigns.loc[i, "Total_consistent"],
           alt_assigns.loc[i, "Sum_prob"]))
    if max_consistent > last_consistent:
        a.make_assign_df(matching_dict[i], set_assign_df=True)
        a.check_assignment_consistency()




# Count correct assignments
(a.assign_df["Res_name"].str.lstrip()==a.assign_df["SS_name"].str.lstrip("_")).sum()

#%% 2nd attempt: more general function for including/excluding potential matches.
consistent = a.assign_df.loc[(a.assign_df["Num_good_links_prev"]==3) & 
                               (a.assign_df["Num_good_links_next"]==3),
                               ["SS_name", "Res_name"]]
inconsistent = a.assign_df.loc[~a.assign_df["Res_name"].isin(consistent["Res_name"]), ["SS_name", "Res_name"]]


b = NAPS_assigner()
b.obs = a.obs.loc[inconsistent["SS_name"],:]
b.preds = a.preds.loc[inconsistent["Res_name"],:]

obs2 = b.obs
preds2 = b.preds

b.add_dummy_rows()
b.calc_log_prob_matrix(sf=1, verbose=False)



def find_matching_constrained(assigner, 
                              inc=pd.DataFrame(columns=["SS_name","Res_name"]), 
                              exc=pd.DataFrame(columns=["SS_name","Res_name"])):
    """Find the best assignment, given constraints
    
    Returns a data frame of SS_names and Res_names of the best matching
    
    inc: a DataFrame of (SS,Res) pairs which must be part of the assignment. 
        First column has the SS_names, second has the Res_names .
    exc: a DataFrame of (SS,Res) pairs which may not be part of the assignment.
    """
    obs = assigner.obs
    preds = assigner.preds
    log_prob_matrix = assigner.log_prob_matrix
    
    # Check that inc and exc lists are consistent
    # Get rid of any entries in exc which chare a Res or SS with inc
    exc_in_inc = exc["SS_name"].isin(inc["SS_name"]) | exc["Res_name"].isin(inc["Res_name"])
    if any(exc_in_inc):
        print("Some values in exc are also found in inc, so are redundant.")
        exc = exc.loc[~exc_in_inc, :]
    # Check for conflicting entries in inc
    conflicts = inc["SS_name"].duplicated(keep=False) | inc["Res_name"].duplicated(keep=False)
    if any(conflicts):
        print("Error: entries in inc conflict with one another")
        print(inc.loc[conflicts,:])
    
    # Removed fixed assignments from probability matrix
    log_prob_matrix_reduced = log_prob_matrix.drop(index=inc["SS_name"]).drop(columns=inc["Res_name"])
    
    # Penalise excluded SS,Res pairs
    penalty = 2*log_prob_matrix.min().min()
    for index, row in exc.iterrows():
        # Need to account for dummy residues or spin systems
        if preds.loc[row["Res_name"], "Dummy_res"]:
            log_prob_matrix_reduced.loc[row["SS_name"], 
                            preds.loc[preds["Dummy_res"],"Res_name"]] = penalty
        elif obs.loc[row["SS_name"], "Dummy_SS"]:
            log_prob_matrix_reduced.loc[obs.loc[obs["Dummy_SS"],"SS_name"], 
                                        row["Res_name"]] = penalty
        else:
            log_prob_matrix_reduced.loc[row["SS_name"], row["Res_name"]] = penalty
    
    row_ind, col_ind = linear_sum_assignment(-1*log_prob_matrix_reduced)
    
    
    # Construct results dataframe
    assignment_reduced = pd.DataFrame({"SS_name":log_prob_matrix_reduced.index[row_ind],
                                       "Res_name":log_prob_matrix_reduced.columns[col_ind]})
#    fixed_assignments = pd.DataFrame({"SS_name":inc, "Res_name":inc_res})

    assignment = pd.concat([inc, assignment_reduced]) 
    return(assignment)

tmp=find_matching_constrained(a)

assign_dict = {}
scores_dict = {}
for i, row in inconsistent.iterrows():
    print(row["Res_name"])
    tmp = find_matching_constrained(a, inc=consistent, exc = inconsistent.loc[[i],:])
    assign_dict[row["SS_name"]] = tmp
    scores_dict[row["SS_name"]] = sum(a.log_prob_matrix.lookup(tmp["SS_name"], tmp["Res_name"]))

#%% 3rd attempt: Rigorously find k-best assignments (don't worry about consistency for now.

from collections import namedtuple
from sortedcontainers import SortedListWithKey

def score_matching(assigner, matching):
    "Calculate the total score of a matching"
    # Calculate sum probability for the best matching
    return(sum(assigner.log_prob_matrix.lookup(matching["SS_name"], 
                                               matching["Res_name"])))
    
Node = namedtuple("Node", ["sum_log_prob","matching","inc","exc"])

# Note: python library says removing items from the start of a list is slow, 
# and that a collections.deque object may be better

matching1.index = matching1["SS_name"]
matching1.index.name = None

ranked_nodes = SortedListWithKey(key=lambda n: n.sum_log_prob)
unranked_nodes = SortedListWithKey([Node(score_matching(a, matching1),
                                         matching1, inc=None, exc=None)],
                                   key=lambda n: n.sum_log_prob)
k=10
while len(ranked_nodes)<k:
    #print(len(ranked_nodes))
    # Set highest scoring unranked node as current_node
    current_node = unranked_nodes.pop()
    
    #Print status
    s = str(len(ranked_nodes))+"\t"
    if current_node.inc is not None:
        s = s + "inc:" +str(len(current_node.inc))+ "\t"
    if current_node.exc is not None:
        s = s + "exc:"+ str(len(current_node.exc))
    print(s)
    # If the current node has forced included pairings, get a list of all non-mandatory pairs
    if current_node.inc is not None:
        matching_reduced = current_node.matching[~current_node.matching["SS_name"].isin(current_node.inc["SS_name"])]
    else:
        matching_reduced = current_node.matching
    
    # Create child nodes and add them to the unranked_nodes list
    for i in range(matching_reduced.shape[0]):
        #print(i)
        #Construct dataframes of excluded and included pairs
        exc_i = matching_reduced.iloc[[i],:]
        if current_node.exc is not None:
            exc_i = exc_i.append(current_node.exc, ignore_index=True)
        inc_i = matching_reduced.iloc[0:i,:]
        if current_node.inc is not None:
            inc_i = inc_i.append(current_node.inc, ignore_index=True)
        inc_i = inc_i.drop_duplicates()
        # If there are no included pairs, it's better for inc_i to be None than an empty df
        if inc_i.shape[0]==0: 
            inc_i = None
        
        matching_i = a.find_best_assignments(inc=inc_i, exc=exc_i, 
                                             return_none_if_all_dummy=True,
                                             verbose=False)
        if matching_i is None:
            pass
        else:
            matching_i.index = matching_i["SS_name"]
            matching_i.index.name = None
            node_i = Node(score_matching(a, matching_i), matching_i, inc_i, exc_i)
            unranked_nodes.add(node_i)
    ranked_nodes.add(current_node)

tmp = [n.matching["Res_name"] for i,n in enumerate(reversed(ranked_nodes))]
results_df = pd.concat(tmp, axis=1)
results_df.columns = list(range(k))
#results_df.index =ranked_nodes[0].matching["SS_name"]
correct_mat = (results_df == results_df.index.tolist()).astype(int)

import matplotlib.pyplot as plt
import seaborn as sns
#plt.pcolormesh(correct_mat.transpose())
#plt.yticks(np.arange(0.5, len(correct_mat.index), 1), correct_mat.index)
#plt.xticks(np.arange(0.5, len(correct_mat.columns), 1), correct_mat.columns)
p = sns.heatmap(correct_mat.transpose())
fig = p.get_figure()
fig.savefig("../plots/test.png", dpi=400)
#plt.show()
#plt.savefig("../plots/test.png")

#tmp = [(n.exc.values[0][0],n.exc.values[0][1]) for n in unranked_nodes]
tmp = [n.sum_log_prob for n in unranked_nodes]
tmp = [n.matching for n in ranked_nodes]

tmp = ranked_nodes[0].matching
tmp2 = ranked_nodes[1].matching
