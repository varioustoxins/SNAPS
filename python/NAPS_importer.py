# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 19:08:55 2018

@author: Alex
"""

import numpy as np
import pandas as pd
import itertools
from pathlib import Path
from Bio.SeqUtils import seq1

class NAPS_importer:
    # Attributes
    peaklists = {}
    assignments = {}
    roots = None
    peaks = None
    obs = None
    
    # Functions
    def import_hsqc_ccpn(self, filename):
        """ Import a CCPN HSQC peaklist.
        
        The peaklist file should have been exported from the CCPN Analysis 
        peak list window.
        If seq_hadamac_details==True, it is assumed the details column contains
        a list of possible amino acid types for the i-1 residue.
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
    
    def import_hsqc_peaks(self, filename, filetype):
        """ Import a peak list for an HSQC spectrum. 
        
        filename: Path to file containing an HSQC peak list.
        filetype: Allowed values are "ccpn", "sparky", "xeasy" or "nmrpipe"
        """
        # Import from file
        if filetype=="ccpn":
            hsqc = pd.read_table(filename) 
            hsqc = hsqc[["Assign F1","Position F1","Position F2", "Height"]]
            hsqc.columns = ["SS_name","H","N","Height"]      
            #hsqc.index = hsqc["SS_name"]
        elif filetype=="sparky":
            hsqc = pd.read_table(filename, sep="\s+")
            hsqc.columns = ["Name","H","N","Height"]
            # If assigned, Name has format "A123HN-A123N"
            # If unassigned, Name column contains "?-?"
            # Get just the first part before the hyphen
            tmp = list(zip(*hsqc["Name"].str.split("-")))[0]
            hsqc["Name"] = tmp
            
            # Deal with assigned peaks first
            # Split atom type from residue name
            tmp = [x.rsplit("H", maxsplit=1) for x in 
                   hsqc.loc[hsqc["Name"]!="?","Name"]]
            # Only keep the backbone amides (Type = "N")
            hsqc.loc[hsqc["Name"]!="?","SS_name"] = list(zip(*tmp))[0]
            hsqc.loc[hsqc["Name"]!="?","Type"] = list(zip(*tmp))[1]
            hsqc.loc[hsqc["Name"]=="?","Type"] = "unknown"
            hsqc = hsqc.loc[hsqc["Type"].isin(["N", "unknown"]),:]
            
            # Now give arbitrary names to the unassigned peaks
            N_unassigned = sum(hsqc["Name"]=="?")
            hsqc.loc[hsqc["Name"]=="?", "SS_name"] = list("x"+
                    (10000+10*pd.Series(range(N_unassigned))).astype(str))
            
            hsqc = hsqc[["SS_name", "N", "H", "Height"]]
        elif filetype=="xeasy":
            hsqc = pd.read_table(filename, sep="\s+", comment="#", header=None,
                                 usecols=[0,1,2,5], 
                                 names=["SS_name","H","N","Height"])
            hsqc["SS_name"] = hsqc["SS_name"].astype(str)
        elif filetype=="nmrpipe":
            # Work out where the column names and data are
            with open(filename, 'r') as f:
                for num, line in enumerate(f, 1):
                    if line.find("VARS")>-1:
                        colnames_line = num
                        colnames = line.split()[1:]
            
            hsqc = pd.read_table(filename, sep="\s+", skiprows=colnames_line+1,
                                names=colnames)
            hsqc = hsqc[["INDEX", "ASS", "X_PPM", "Y_PPM", "HEIGHT"]]
            hsqc.columns = ["ID", "SS_name", "H", "N", "Height"]
            # Use ASS as SS_name if available, otherwise use ID
            hsqc.loc[hsqc["SS_name"]!=None,"SS_name"] = \
                hsqc.loc[hsqc["SS_name"]!=None,"ID"].astype(str)
            hsqc = hsqc[["SS_name", "H", "N", "Height"]]
        else:
            print("import_hsqc_peaks: invalid filetype '%s'." % (filetype))
            return(None)
        
        # Check whether columns go mixed up somehow, and correct if so
        if hsqc["H"].mean() > hsqc["N"].mean():
            tmp = hsqc["N"]
            hsqc["N"] = hsqc["H"]
            hsqc["H"] = tmp
        
        hsqc.index = hsqc["SS_name"]
        
        self.peaklists["hsqc"] = hsqc
        self.roots = hsqc
        
        return(self.roots)
    
    def import_3d_peaks_ccpn(self, filename, spectrum):
        """ Import a CCPN 3D peak list
        
        The peak list should have been exported from the CCPN Analysis peak 
        list dialogue.
        spectrum can be one of "hnco", "hncaco", "hnca", "hncoca", "hncacb" or
        "hncocacb"
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
        
        # Check that the correct dimensions were identified
        # Also check that spectrum argument is valid
        if spectrum in ["hnco", "hncaco", "hnca","hncoca","hncacb","hncocacb"]:
            if set(dim.keys()) != set(["H","N","C"]):
                print("Error: couldn't identify "+spectrum+" columns.") 
        else:
            print("Invalid value of argument: spectrum.")
        
        # Choose and rename columns. 
        # Also, only keep spin systems that are in self.roots
        peaks = peaks[["Assign "+dim["H"], "Position "+dim["H"], 
                       "Position "+dim["N"], "Position "+dim["C"], "Height"]]
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
            
    def import_obs_shifts(self, filename, filetype, SS_num=False):
        """ Import a chemical shift list
        
        filename: Path to text file containing chemical shifts.
        filetype: Allowed values are "ccpn", "sparky", "xeasy" or "nmrpipe"
            (The ccpn option is for importing a Resonance table exported from 
            Analysis v2.x)
        SS_num: If true, will extract the longest number from the SS_name and 
        treat it as the residue number. Without this, it is not possible to get
        the i-1 shifts for each spin system.
            
        """
        # Import from file
        if filetype=="ccpn":
            obs = pd.read_table(filename)
            obs = obs.loc[:,["Residue", "Assign Name", "Shift"]]
            obs.columns = ["SS_name", "Atom_type", "Shift"]
            obs["Atom_type"] = obs["Atom_type"].str.upper()
        elif filetype=="sparky":
            obs = pd.read_table(filename, sep="\s+")
            obs = obs.loc[:,["Group", "Atom", "Shift"]]
            obs.columns = ["SS_name", "Atom_type", "Shift"]
            obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
        elif filetype=="xeasy":
            obs = pd.read_table(filename, sep="\s+", 
                                header=None, na_values="999.000",
                                names=["i","Shift","SD","Atom_type","SS_name"])
            obs = obs.loc[:, ["SS_name", "Atom_type", "Shift"]]
            obs["SS_name"] = obs["SS_name"].astype(str)
            obs = obs.dropna(subset=["Shift"])
            obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
        elif filetype=="nmrpipe":
            # Work out where the column names and data are
            with open(filename, 'r') as f:
                for num, line in enumerate(f, 1):
                    if line.find("VARS")>-1:
                        colnames_line = num
            
            obs = pd.read_table(filename, sep="\s+", skiprows=colnames_line+1, 
                                names=["SS_name","Res_type","Atom_type","Shift"])
            obs = obs.loc[:, ["SS_name", "Atom_type", "Shift"]]
            obs["SS_name"] = obs["SS_name"].astype(str)
            obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
        else:
            print("import_obs_shifts: invalid filetype '%s'." % (filetype))
            return(None)
        
        # Restrict to backbone atom types
        obs = obs.loc[obs["Atom_type"].isin(["H","HA","N","C","CA","CB",
                                             "Cm1","CAm1","CBm1"]),:]
        
        # Convert from long to wide
        obs = obs.pivot(index="SS_name", columns="Atom_type", values="Shift")
        obs.insert(0, "SS_name", obs.index.values)
        
        # Extract residue number from SS_name (if present), and get m1 shifts
        if SS_num:
            obs["Res_N"] = obs["SS_name"].str.extract(r"(\d+)").astype(int)
            
            obs.index = obs["Res_N"]
            obs_m1 = obs[list({"C","CA","CB"}.intersection(obs.columns))]
            obs_m1.index = obs_m1.index+1
            obs_m1.columns = obs_m1.columns + "m1"
            obs = pd.merge(obs, obs_m1, how="left", 
                           left_index=True, right_index=True)
            obs = obs.drop(columns="Res_N")
        
        self.obs = obs
        return(self.obs)
        
        
    
                    

#%%


def find_all_assignments(peaks, atoms):
    """ Find all possible assignments of peaks to atoms, including where some or all atoms are unassigned
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

#path = "/Users/aph516/GitHub/NAPS/data/P3a_L273/"
path = "C:/Users/Alex/GitHub/NAPS/data/P3a_L273R/"

a = NAPS_importer()

a.import_hsqc_ccpn(path+"hsqc HADAMAC.txt")
a.import_3d_peaks_ccpn(path+"hnco.txt", "hnco")
a.import_3d_peaks_ccpn(path+"cbcaconh.txt", "hncocacb")
a.import_3d_peaks_ccpn(path+"hncacb.txt", "hncacb")
#a.guess_peak_atom_types()

tmp = a.import_hsqc_peaks(path+"nmrpipe_peaks.tab", "nmrpipe")

tmp = a.import_obs_shifts(path+"ccpn_resonances.txt", "ccpn", SS_num=True)
tmp = a.import_obs_shifts(path+"sparky shifts.txt", "sparky", SS_num=True)
tmp = a.import_obs_shifts(path+"Xeasy shifts.txt", "xeasy", SS_num=True)
tmp = a.import_obs_shifts(path+"nmrpipe_shifts.tab", "nmrpipe", SS_num=True)

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
