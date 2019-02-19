#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:30:36 2018

@author: aph516
"""

#import pandas as pd
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
import argparse
#from pathlib import Path
import logging

#### User input

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
#parser.add_argument("--delta_correlation", action="store_true", 
#                    help="If set, account for correlations between prediction errors of different atom types")
parser.add_argument("-a", "--alt_assignments", default=-1, type=int,
                    help="The number of alternative assignments to generate, "+
                    "in addition to the highest ranked.")
parser.add_argument("--plot_stem", 
                    default="/Users/aph516/GitHub/NAPS/plots/plot",
                    help="The base name for plots (do not include a file extension).")

if True:
    args = parser.parse_args()
else:
    args = parser.parse_args(("../data/P3a_L273R/naps_shifts.txt",
                              "../output/test.txt",
                              "--shift_type","naps",
                              "--pred_file","../data/P3a_L273R/shiftx2.cs",
                              "-c","../config/config.txt",
                              "-l","../output/test.log"))    

# Set up logging
if isinstance(args.log_file, str):
    print(args.log_file)
    logging.basicConfig(filename=args.log_file, filemode='w', level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s")
else:
    logging.basicConfig(level=logging.ERROR)

#%%
a = NAPS_assigner()

# Import config file
a.read_config_file(args.config_file)
logging.info("Read in configuration from %s.", args.config_file)

# Account for any command line arguments that overide config file
if args.alt_assignments>=0:
    a.pars["alt_assignments"] = args.alt_assignments
#if args.delta_correlation:
#    a.pars["prob_method"] = "delta_correlation"

# Import observed and predicted shifts
importer = NAPS_importer()

if args.shift_type=="test":
    importer.import_testset_shifts(args.input_file)
else:
    importer.import_obs_shifts(args.input_file, args.shift_type, SS_num=False)
a.obs = importer.obs
logging.info("Read in %d spin systems from %s.", 
             len(a.obs["SS_name"]), args.input_file)


a.import_pred_shifts(args.pred_file, args.pred_type)
logging.info("Read in %d predicted residues from %s.", 
             len(a.preds["Res_name"]), args.pred_file)

# Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix2(sf=1, verbose=False)
logging.info("Calculated log probability matrix (%dx%d).", 
             a.log_prob_matrix.shape[0], a.log_prob_matrix.shape[1])
matching = a.find_best_assignments()
a.make_assign_df(matching, set_assign_df=True)
logging.info("Calculated best assignment.")
assign_df = a.check_assignment_consistency(threshold=0.1)
logging.info("Checked assignment consistency.")

if a.pars["alt_assignments"]>0:
    a.find_alt_assignments(N=a.pars["alt_assignments"], verbose=False, 
                            by_ss=True)
    logging.info("Calculated the %d next best assignments for each spin system", 
                 a.pars["alt_assignments"])
    a.alt_assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                           index=False)
else:
    a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f", 
                       index=False)
    
logging.info("Wrote results to %s", args.output_file)
#tmp = a.find_alt_assignments(best_match_indexes, by_res=False)
#tmp = a.find_alt_assignments2(N=2, verbose=True, by_ss=False)

# Make some plots
if a.pars["plot_strips"]:
    plt = a.plot_strips()
    plt.save(args.plot_stem+"_strips.pdf", height=210, width=max(297,297/80*a.assign_df["SS_name"].count()), units="mm")
    logging.info("Wrote strip plot to %s", args.plot_stem+"_strips.pdf")