#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 10:45:35 2018

@author: aph516
"""

import numpy as np
import pandas as pd

# Import the observed chemical shifts
obs = pd.read_table("~/GitHub/NAPS/data/testset/simplified_BMRB/4032.txt")

#### Import the predicted chemical shifts
preds = pd.read_csv("~/GitHub/NAPS/data/testset/shiftx2_results/A001_1KF3A.cs")
preds["Res_name"] = preds["NUM"].astype(str)+preds["RES"]
preds = preds.reindex(columns=["NUM","RES","Res_name","ATOMNAME","SHIFT"])
preds.columns = ["Res_N","Res_type","Res_name","Atom_type","Shift"]

# Convert from wide to long format
preds_wide = preds.pivot(index="Res_N", columns="Atom_type", values="Shift")

# Make columns for the i-1 predicted shifts of C, CA and CB
preds_wide_m1 = preds_wide[["C","CA","CB"]].copy()
preds_wide_m1.index = preds_wide_m1.index+1
#preds_wide_m1 = preds_wide_m1[["C","CA","CB"]]
preds_wide_m1.columns = ["Cm1","CAm1","CBm1"]
preds_wide = pd.concat([preds_wide, preds_wide_m1], axis=1)


# Restrict to only certain atom types
atom_list = ["H","N","C","CA","CB","Cm1","CAm1","CBm1","HA"]
preds_wide = preds_wide[atom_list]
preds_wide