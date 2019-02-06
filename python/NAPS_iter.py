# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 21:32:56 2019

@author: Alex
"""

import pandas as pd
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
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
config = pd.read_table(args.config_file, sep="\s+", comment="#", header=None,
                       index_col=0, names=["Value"])
a.pars["pred_offset"] = int(config.loc["pred_offset"].Value)
a.pars["prob_method"] = config.loc["prob_method"].Value
a.pars["alt_assignments"] = int(config.loc["alt_assignments"].Value)
a.pars["atom_set"] = {s.strip() for s in config.loc["atom_set"].Value.split(",")}
tmp = [s.strip() for s in config.loc["atom_sd"].Value.split(",")]
a.pars["atom_sd"] = dict([(x.split(":")[0], float(x.split(":")[1])) for x in tmp])
a.pars["plot_strips"] = bool(strtobool(config.loc["plot_strips"].Value))
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
a.calc_log_prob_matrix(sf=1, verbose=False)
logging.info("Calculated log probability matrix (%dx%d).", 
             a.log_prob_matrix.shape[0], a.log_prob_matrix.shape[1])
assign_df, best_match_indexes = a.find_best_assignment()
logging.info("Calculated best assignment.")
assign_df = a.check_assignment_consistency(threshold=0.1)
logging.info("Checked assignment consistency.")

obs = a.obs
preds = a.preds

# Work out which residues are consistent
consistent_res = assign_df.loc[(assign_df["Num_good_links_prev"]==3) & (assign_df["Num_good_links_next"]==3),"Res_name"]
inconsistent_res = assign_df.loc[~assign_df["Res_name"].isin(consistent_res), "Res_name"]
assign_df.index = assign_df["Res_name"]
inconsistent_SS = assign_df.loc[inconsistent_res, "SS_name"]

b = NAPS_assigner()
b.obs = a.obs.loc[inconsistent_SS,:]
b.preds = a.preds.loc[inconsistent_res,:]

obs2 = b.obs
preds2 = b.preds

b.add_dummy_rows()
b.calc_log_prob_matrix(sf=1, verbose=False)

def find_matching_constrained(assigner, inc=None, exc=None):
    """Find the best assignment, given constraints
    
    Returns a data frame of SS_names and Res_names of the best matching
    
    inc: a list of (SS, Res) tuples which must be part of the assignment.
    exc: a list of (SS, Res) tuples which may not be part of the assignment.
    """
    obs = assigner.obs
    preds = assigner.preds
    log_prob_matrix = assigner.log_prob_matrix
    
    # Check that inc and exc lists are consistent
    
    # Removed fixed assignments from probability matrix
    inc_SS, inc_res = zip(*inc)     # Make lists of SS and residues that are fixed
    log_prob_matrix_reduced = log_prob_matrix.drop(index=list(inc_SS)).drop(columns=list(inc_res))
    
    # Penalise excluded SS,Res pairs
    
    # Construct results dataframe
    fixed_assignments = pd.DataFrame({"SS_name":inc_SS, "Res_name":inc_res})    
    pass

#assign_df, best_match_indexes = b.find_best_assignment()
#assign_df = b.check_assignment_consistency(threshold=0.1)

# Make a strip plot
plt = a.plot_strips()
plt.save(args.plot_file, height=210, width=max(297,297/80*a.assign_df["SS_name"].count()), units="mm", limitsize=False)
logging.info("Wrote strip plot to %s", args.plot_file)

#if a.pars["alt_assignments"]>0:
#    a.find_alt_assignments2(N=a.pars["alt_assignments"], verbose=False, 
#                            by_ss=True)
#    logging.info("Calculated the %d next best assignments for each spin system", 
#                 a.pars["alt_assignments"])
#    a.alt_assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
#                           index=False)
#    logging.info("Wrote results to %s", args.output_file)
#else:
#    a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
#                       index=False)
#    logging.info("Wrote results to %s", args.output_file)


