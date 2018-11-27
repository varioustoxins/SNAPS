#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:30:36 2018

@author: aph516
"""

import pandas as pd
from NAPS_assigner import NAPS_assigner
import argparse

#### User input

# For testing
if False:
    runfile('/Users/aph516/GitHub/NAPS/python/NAPS.py', wdir='/Users/aph516/GitHub/NAPS/python',
            args="../data/testset/simplified_BMRB/4032.txt "+
            "../data/testset/shiftx2_results/A001_1KF3A.cs "+
            "../output/test.txt")
    runfile('/Users/aph516/GitHub/NAPS/python/NAPS.py', wdir='/Users/aph516/GitHub/NAPS/python',
            args="../data/testset/simplified_BMRB/4032.txt "+
            "../data/testset/shiftx2_results/A001_1KF3A.cs "+
            "../output/test.txt"+
            "-c config.txt")

parser = argparse.ArgumentParser(description="NMR Assignments from Predicted Shifts")
parser.add_argument("shift_file")
parser.add_argument("pred_file")
parser.add_argument("out_file")
parser.add_argument("-c", "--config_file", default="/Users/aph516/GitHub/NAPS/python/config.txt")


args = parser.parse_args()

print(args.config_file)




#%%



#### Set up the assigner

a = NAPS_assigner()

# Import config file
config = pd.read_table("/Users/aph516/GitHub/NAPS/python/config.txt", sep="\s+", comment="#", header=None,
                       index_col=0, names=["Value"])
a.pars["shiftx2_offset"] = int(config.loc["shiftx2_offset"].Value)
a.pars["atom_set"] = {s.strip() for s in config.loc["atom_set"].Value.split(",")}
tmp = [s.strip() for s in config.loc["atom_sd"].Value.split(",")]
a.pars["atom_sd"] = dict([(x.split(":")[0], float(x.split(":")[1])) for x in tmp])


a.import_obs_shifts(args.shift_file)
a.read_shiftx2(args.pred_file)

#### Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix(sf=2, verbose=False)
assign_df, best_match_indexes = a.find_best_assignment()
a.check_assignment_consistency(threshold=0.1)
a.assign_df.to_csv(args.out_file, sep="\t", float_format="%.3f")

#tmp = a.find_alt_assignments(best_match_indexes, by_res=False)
tmp = a.find_alt_assignments2(N=2, verbose=True)


#%%

#### Make some plots