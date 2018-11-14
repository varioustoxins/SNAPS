#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 11:00:12 2018
Test new method for getting shifts from an unassigned peak list.
Try all possible assignments, elliminate impossible ones, then rank by plausibility.

@author: aph516
"""

import numpy as np
import pandas as pd
import itertools

peaks = [1,2,3]
atoms = ['a', 'b', 'c']


a = list(itertools.permutations(peaks, 2))  # Use to generate all the possible orderings
b = list(itertools.combinations(peaks, 2))  # Could use to choose which atoms are missing

# Start assuming no atoms are missing
df0 = pd.DataFrame(list(itertools.permutations(peaks, 3)), columns=atoms)

# Now with 1 missing atom
# missing = list(itertools.combinations(peaks, 1)) # No point for one missing
df1 = pd.DataFrame(columns=atoms)
tmp = pd.DataFrame(list(itertools.permutations(peaks, 2)))
for a in atoms:
#    tmp2 = tmp
#    tmp2.columns = [x for x in atoms if x is not a]
#    tmp2[a] = np.NaN
    tmp2 = tmp.copy()
    tmp2.columns = [x for x in atoms if x is not a]
    df1 = df1.append(tmp2)
    #pd.concat(df1, tmp2)

# Now with 2 missing atoms
df2 = pd.DataFrame(columns=atoms)
tmp = pd.DataFrame(list(itertools.permutations(peaks, 1)))
missing = list(itertools.combinations(atoms, 2))
for m in missing:
    tmp2 = tmp.copy()
    tmp2.columns = [x for x in atoms if x not in m]
    df2 = df2.append(tmp2)
    
# With 3 missing
df3 = pd.DataFrame(columns=atoms)
df3.loc[0] = [np.NaN]*3

df = pd.concat([df0, df1, df2, df3], ignore_index=True)


# Repeat for arbitrary numpers of peaks and atoms
def find_all_assignments(peaks, atoms):
    """
    Find all possible assignments of peaks to atoms, including where some or all atoms are unassigned
    """
    if len(atoms)==0:   # Case where atoms list is empty
        return(None)
        
    df = pd.DataFrame(columns=atoms)
    
    if len(peaks)==0:   # Case where peaks list is empty
        return(df)
    
    # Generate all possible assignments for every possible u (the number of unassigned atoms)
    u_min = max(0, len(atoms) - len(peaks))
    for u in range(u_min, len(atoms)+1):
        df2 = pd.DataFrame(columns=atoms)
        tmp = pd.DataFrame(list(itertools.permutations(peaks, len(atoms)-u)))
        missing = list(itertools.combinations(atoms, u))
        for m in missing:
            tmp2 = tmp.copy()
            tmp2.columns = [x for x in atoms if x not in m]
            df2 = df2.append(tmp2)
        df = df.append(df2, ignore_index=True)
    return(df)
    
a = find_all_assignments([], [])
b = find_all_assignments([], ['A','B'])
c = find_all_assignments([1,2], [])
d = find_all_assignments([1], ['A','B'])
e = find_all_assignments([1,2], ['A'])
f = find_all_assignments([1,2,3], ['A','B','C'])
