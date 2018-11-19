# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 19:08:55 2018

@author: Alex
"""

import numpy as np
import pandas as pd
import itertools
from Bio.SeqUtils import seq1

class NAPS_importer:
    # Attributes
    peaklists = {}
    assignments = {}
    roots = None
    peaks = None
    shifts = None
    
    # Functions
    def import_hsqc_ccpn(self, filename):
        """ Import an HSQC peaklist, as exported from the CCPN Analysis peak list dialogue.
        
        If seq_hadamac_details==True, it is assumed the details column contains a list of possible amino acid types for the i-1 residue.
        """
        # TODO: Check that SS names are unique
        
        hsqc = pd.read_table(filename) 
        
        # Keep just the columns we need, and rename
        hsqc = hsqc[["Assign F1","Position F1","Position F2", "Height"]]
        hsqc.columns = ["SS_name","H","N","Height"]      
        hsqc.index = hsqc["SS_name"]
        
        self.peaklists["hsqc"] = hsqc
        self.roots = hsqc
        
        return (self)
    
    def import_3d_peaks_ccpn(self, filename, spectrum):
        """ Import a 3D peaklist, as exported from the CCPN Analysis peak list dialogue.
        
        spectrum can be one of "hnco", "hncaco", "hnca", "hncoca", "hncacb", "hncocacb"
        """
        peaks = pd.read_table(filename)
        
        # Work out which column contains which dimension
        dim = {}
        for f in ["F1","F2","F3"]:
            if peaks["Position "+f].mean() > 150:
                dim["C"] = f
            elif peaks["Position "+f].mean() > 100:
                dim["N"] = f
            elif peaks["Position "+f].mean() > 15:
                dim["C"] = f
            elif peaks["Position "+f].mean() > 6:
                dim["H"] = f
            elif peaks["Position "+f].mean() > 0:
                dim["HA"] = f
            else: 
                dim["Unknown"] = f
        
        # Check that the correct dimensions were identified, and also that spectrum argument is valid
        if spectrum in ["hnco", "hncaco", "hnca","hncoca","hncacb","hncocacb"]:
            if set(dim.keys()) != set(["H","N","C"]):
                print("Error: couldn't identify "+spectrum+" columns.") 
        else:
            print("Invalid value of argument: spectrum.")
        
        # Choose and rename columns. Also, only keep spin systems that are in self.roots
        peaks = peaks[["Assign "+dim["H"], "Position "+dim["H"], "Position "+dim["N"], "Position "+dim["C"], "Height"]]
        peaks.columns = ["SS_name", "H", "N", "C", "Height"]
        peaks["Spectrum"] = spectrum
        peaks = peaks.loc[peaks["SS_name"].isin(self.roots["SS_name"])]
        
        # Add peaks to the peaklist dataframe
        self.peaklists[spectrum] = peaks
#        if isinstance(self.peaks, pd.DataFrame):
#            self.peaks = self.peaks.append(peaks, ignore_index=True)
#        else:
#            self.peaks = peaks
            
        return(self)
    
    def guess_peak_atom_types(self):
        """ Makes a best guess of the atom type for all peaks in self.peaks
        """
        
        self.peaks["Atom_type"] = "NA"
        
        for ss in self.roots["SS_name"]:
            ss_peaks = self.peaks.loc[self.peaks["SS_name"]==ss,:]
            
            if "hnco" in ss_peaks["Spectrum"].values:
                # Take the strongest HNCO peak, and set its atom type to Cm1
                i = ss_peaks.loc[ss_peaks["Spectrum"]=="hnco","Height"].idxmax()
                ss_peaks.loc[i,"Atom_type"] = "Cm1"

            self.peaks.loc[self.peaks["SS_name"]==ss,:] = ss_peaks   
                    

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

def score_plausibility(assignments):
    atom_list = {"CA","CB","CAm1","CBm1"}.intersection(assignments.columns)
    scores = pd.DataFrame(index = assignments.index, columns=["Rest","Gly","Ser","Thr"])
    scores.iloc[:,:] = 0
    tmp = pd.DataFrame(index = assignments.index)
    for a in atom_list:
        if a in ["CA", "CAm1"]:
            scores.loc[assignments[a].between(48,68),"Rest"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(46,70) & assignments[a].notnull(),"Rest"] = -2
            
            scores.loc[assignments[a].between(43,48),"Gly"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(41,50) & assignments[a].notnull(),"Gly"] = -2
            
            scores.loc[assignments[a].between(54,63),"Ser"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(51,67) & assignments[a].notnull(),"Ser"] = -2
            
            scores.loc[assignments[a].between(57,69),"Thr"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(52,72) & assignments[a].notnull(),"Thr"] = -2
            
        elif a in ["CB", "CBm1"]:
            scores.loc[assignments[a].between(15,45),"Rest"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(13,53) & assignments[a].notnull(),"Rest"] = -2
            
            scores.loc[assignments[a].between(999,999),"Gly"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(999,999) & assignments[a].notnull(),"Gly"] = -2
            
            scores.loc[assignments[a].between(71,67),"Ser"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(58,69) & assignments[a].notnull(),"Ser"] = -2
            
            scores.loc[assignments[a].between(66,73),"Thr"] += 1/len(atom_list)
            scores.loc[~assignments[a].between(64,76) & assignments[a].notnull(),"Thr"] = -2

    return(scores)

#%%

#### Test code

a = NAPS_importer()
a.import_hsqc_ccpn("~/GitHub/NAPS/data/P3a_L273R/hsqc HADAMAC.txt")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/hnco.txt", "hnco")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/cbcaconh.txt", "hncocacb")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/hncacb.txt", "hncacb")
#a.guess_peak_atom_types()

roots = a.roots
peaklists = a.peaklists

hncocacb = a.peaklists["hncocacb"]

# Make a list of all possible peak assignments for the HNCOCACB
atoms = ["CAm1","CBm1"]
assignments = pd.DataFrame(columns=["SS_name","As_ID"]+atoms)
for SS in a.roots["SS_name"]:
    peaks = a.peaklists["hncocacb"].loc[a.peaklists["hncocacb"]["SS_name"]==SS, "C"]
    tmp = find_all_assignments(peaks, atoms)
    tmp["SS_name"] = SS
    tmp["As_ID"] = range(tmp.shape[0])
    assignments = assignments.append(tmp, ignore_index=True)

# Check which assignments are plausible
b = pd.concat([assignments, score_plausibility(assignments)], axis=1)

c = b.loc[b[["Rest","Gly","Ser","Thr"]].apply(max, axis=1)>=0,:]
