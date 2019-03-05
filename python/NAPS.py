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

#### Command line arguments
parser = argparse.ArgumentParser(description="NAPS (NMR Assignments from Predicted Shifts)")
parser.add_argument("shift_file", 
                    help="A table of observed chemical shifts.")
parser.add_argument("pred_file",
                    help="A table of predicted chemical shifts.")
parser.add_argument("output_file",
                    help="The file results will be written to.")

parser.add_argument("--shift_type", 
                    choices=["naps", "ccpn", "sparky", 
                             "xeasy", "nmrpipe", "test"], 
                    default="naps", 
                    help="The format of the observed shift file.")
parser.add_argument("--pred_type", 
                    choices=["shiftx2", "sparta+"],
                    default="shiftx2", 
                    help="The format of the predicted shifts")
parser.add_argument("-c", "--config_file", 
                    default="/Users/aph516/GitHub/NAPS/python/config.txt",
                    help="A file containing parameters for the analysis.")
parser.add_argument("-l", "--log_file", default=None,
                    help="A file logging information will be written to.")
#parser.add_argument("--delta_correlation", action="store_true", 
#                    help="If set, account for correlations between prediction errors of different atom types")
parser.add_argument("-a", "--alt_assignments", default=-1, type=int,
                    help="The number of alternative assignments to generate, "+
                    "in addition to the highest ranked.")
parser.add_argument("--test_aa_classes", default=None, 
                    help="""For test data only. 
                    A string containing a comma-separated list of the amino acid 
                    classes for the i residue, a semicolon, then a list of AA 
                    classes for the i-1 residue. No spaces. 
                    eg. "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML" for 
                    a sequential HADAMAC """)
parser.add_argument("--plot_file", 
                    default="/Users/aph516/GitHub/NAPS/plots/plot",
                    help="A filename for any output plots.")

if True:
    args = parser.parse_args()
else:   # For testing
    args = parser.parse_args(("../data/P3a_L273R/naps_shifts.txt",
                              "../output/test.txt",
                              "--shift_type","naps",
                              "--pred_file","../data/P3a_L273R/shiftx2.cs",
                              "-c","../config/config.txt",
                              "-l","../output/test.log"))    

# Set up logging
if isinstance(args.log_file, str):
    logging.basicConfig(filename=args.log_file, filemode='w', level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s")
else:
    logging.basicConfig(level=logging.ERROR)

#%%
#### Set up the NAPS_assigner object
a = NAPS_assigner()

# Import config file
a.read_config_file(args.config_file)
logging.info("Read in configuration from %s.", args.config_file)

# Account for any command line arguments that overide config file
if args.alt_assignments>=0:
    a.pars["alt_assignments"] = args.alt_assignments

# Import observed and predicted shifts
importer = NAPS_importer()

if args.shift_type=="test":
    if args.test_aa_classes is None:
        importer.import_testset_shifts(args.shift_file)
    else:
        AA_class, AA_classm1 = args.test_aa_classes.split(";")
        importer.import_testset_shifts(args.shift_file,
                                       SS_class=AA_class.split(","),
                                       SS_classm1=AA_classm1.split(","))
else:
    importer.import_obs_shifts(args.shift_file, args.shift_type, SS_num=False)
a.obs = importer.obs
logging.info("Read in %d spin systems from %s.", 
             len(a.obs["SS_name"]), args.shift_file)


a.import_pred_shifts(args.pred_file, args.pred_type)
logging.info("Read in %d predicted residues from %s.", 
             len(a.preds["Res_name"]), args.pred_file)

#### Do the analysis
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

#### Make some plots
if a.pars["plot_strips"]:
    plt = a.plot_strips_bokeh(args.plot_file, "html")
#    plt.save(args.plot_file, height=210, width=max(297,297/80*a.assign_df["SS_name"].count()), 
#             units="mm", limitsize=False)
    logging.info("Wrote strip plot to %s", args.plot_file)