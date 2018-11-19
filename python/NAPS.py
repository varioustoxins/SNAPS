#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:30:36 2018

@author: aph516
"""

from NAPS_assigner import NAPS_assigner
import argparse

#### User input

parser = argparse.ArgumentParser(description="NMR Assignments from Predicted Shifts")
parser.add_argument("shift_file")
parser.add_argument("pred_file")
parser.add_argument("out_file")


args = parser.parse_args()

print(args.shift_file)
print(args.pred_file)
print(args.out_file)




#%%

#### Do the analysis

#a = NAPS_assigner()
#a.import_obs_shifts(args.shift_file)
#a.read_shiftx2(args.pred_file, offset=208)
#a.add_dummy_rows()
#a.calc_log_prob_matrix(sf=2, verbose=False)
#assign_df, best_match_indexes = a.find_best_assignment()



#%%

#### Make some plots