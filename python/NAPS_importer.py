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
from math import sqrt
#import nmrstarlib

class NAPS_importer:
    # Attributes
#    peaklists = {}
#    assignments = {}
#    roots = None
#    peaks = None
#    obs = None
    
    def __init__(self):
        self.peaklists = {}
        #assignments = {}
        self.roots = None
        self.obs = None
    
    # Functions
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
            hsqc.columns = ["SS_name","H","N","Height"]
            # If assigned, Name has format "A123HN-A123N"
            # If unassigned, Name column contains "?-?"
            # Get just the first part before the hyphen
            tmp = list(zip(*hsqc["SS_name"].str.split("-")))[0]
            hsqc["SS_name"] = tmp
            
            # Give arbitrary names to the unassigned peaks
            N_unassigned = sum(hsqc["SS_name"]=="?")
            hsqc.loc[hsqc["SS_name"]=="?", "SS_name"] = list("x"+
                    pd.Series(range(N_unassigned)).astype(str))
            
        elif filetype=="xeasy":
            hsqc = pd.read_table(filename, sep="\s+", comment="#", header=None,
                                 usecols=[0,1,2,5], 
                                 names=["SS_name","H","N","Height"])
            hsqc["SS_name"] = "x" + hsqc["SS_name"].astype(str)
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
        
    def import_3d_peaks(self, filename, filetype, spectrum, 
                        assign_nearest_root=False):
        """Import a 3D peak list in various formats
        
        filetype: one of "ccpn", "sparky", "xeasy" or "nmrpipe"
        spectrum: one of "hnco", "hncaco", "hnca", "hncoca", "hncacb",
            "hncocacb" or "hnha"
        assign_nearest_root: If True, this will assign each peak to the spin 
            system with the closest root (hsqc) peak. If False, peaks will be 
            assigned to spin systems based on information in the original file. 
            In this case, the proton assignment is used for CCPN and Sparky, 
            while the ASS column is used for nmrPipe. Xeasy peaklists alone 
            do not seem to contain assignment information.
        """
        if filetype == "ccpn":
            peaks = pd.read_table(filename,
                                  usecols=["Position F1","Position F2",
                                           "Position F3","Assign F1",
                                           "Assign F2","Assign F3",
                                           "Height"])
            peaks.columns = ["F1","F2","F3","A1","A2","A3","Height"]
        elif filetype == "sparky":
            peaks = pd.read_table(filename, sep="\s+")
            peaks.columns=["Name","F1","F2","F3","Height"]
            peaks["A1"], peaks["A2"], peaks["A3"] = list(zip(
                                                *peaks["Name"].str.split("-")))
            #return(peaks)
        elif filetype == "xeasy":
            peaks = pd.read_table(filename, sep="\s+", comment="#", header=None,
                                 usecols=[1,2,3,6], 
                                 names=["F1","F2","F3","Height"])
            peaks["SS_name"] = None
        elif filetype == "nmrpipe":
            # Work out where the column names and data are
            with open(filename, 'r') as f:
                for num, line in enumerate(f, 1):
                    if line.find("VARS")>-1:
                        colnames_line = num
                        colnames = line.split()[1:]
                        
            peaks = pd.read_table(filename, sep="\s+", skiprows=colnames_line+1,
                                names=colnames)
            peaks = peaks[["ASS", "X_PPM", "Y_PPM", "Z_PPM", "HEIGHT"]]
            peaks.columns = ["SS_name", "F1", "F2", "F3", "Height"]
            
        else:
            print("import_3d_peaks: invalid filetype '%s'." % (filetype))
            return(None)
            
        # Work out which column contains which dimension
        dim = {}
        for i in ["1","2","3"]:
            if peaks["F"+i].mean() > 150:
                dim["C"] = i
            elif peaks["F"+i].mean() > 100:
                dim["N"] = i
            elif peaks["F"+i].mean() > 15:
                dim["C"] = i
            elif peaks["F"+i].mean() > 6:
                dim["H"] = i
            elif peaks["F"+i].mean() > 0:
                dim["HA"] = i
            else: 
                dim["Unknown"] = i
        
        # Check that the correct dimensions were identified
        # Also check that spectrum argument is valid
        if spectrum in ["hnco", "hncaco", "hnca","hncoca","hncacb","hncocacb"]:
            if set(dim.keys()) != set(["H","N","C"]):
                print("Error: couldn't identify "+spectrum+" columns.") 
        elif spectrum == "hnha":
            if set(dim.keys()) != set(["H","N","HA"]):
                print("Error: couldn't identify "+spectrum+" columns.")
        else:
            print("Invalid value of argument: spectrum.")
        
        # Populate SS_name column
        if filetype in ["ccpn", "sparky"]:
            peaks["SS_name"] = peaks["A"+dim["H"]]
        elif filetype == "nmrpipe":
            peaks["SS_name"] = peaks["ASS"]

        # Choose and rename columns. 
        peaks = peaks[["SS_name", "F"+dim["H"], "F"+dim["N"], "F"+dim["C"], 
                       "Height"]]
        peaks.columns = ["SS_name", "H", "N", "C", "Height"]
        peaks["Spectrum"] = spectrum
        
        # If assign_nearest_root, find closest root resonance for each peak 
        # and set that as SS_name.
        if assign_nearest_root or spectrum=="xeasy":
            peaks["SS_name"] == None
            
            roots = self.roots
            for i in peaks.index:
                delta = ((roots["H"]-peaks.loc[i,"H"])**2 
                         + (0.2*(roots["N"]-peaks.loc[i,"N"]))**2).apply(sqrt)
                peaks.loc[i, "SS_name"] = roots.loc[delta.idxmin(), "SS_name"]

        # Also, only keep spin systems that are in self.roots
        #peaks = peaks.loc[peaks["SS_name"].isin(self.roots["SS_name"])]
        
        # Add peaks to the peaklist dictionary
        self.peaklists[spectrum] = peaks
            
        return(self.peaklists[spectrum])
        
    def find_shifts_from_peaks(self):
        """ Work out chemical shifts for each spin system from peak lists
        
        Will use all spectra in the self.peaklists dictionary. The following
        spectra are supported: hnco, hncaco, hnca, hncoca, hncacb, hncocacb.
        """
        # Possible extension: could put a parameter to control hncacb sign interpretation
        obs = pd.DataFrame({"SS_name": self.roots["SS_name"], 
                            "H": self.roots["H"],
                            "N": self.roots["N"]})
        obs.index = obs["SS_name"]
        for spec in self.peaklists.keys():
            peaks = self.peaklists[spec]
            for ss in obs["SS_name"]:
                if ss in peaks["SS_name"].values:
                    ss_peaks = peaks.loc[peaks["SS_name"]==ss,:]
                    if spec=="hnco":
                        # Set Cm1 shift to strongest peak in this spin system.         
                        i = ss_peaks.loc[:,"Height"].idxmax()
                        obs.loc[ss, "Cm1"] = ss_peaks.loc[i, "C"]
                    elif spec=="hncaco":
                        # Set C shift to strongest peak in this spin system.         
                        i = ss_peaks.loc[:,"Height"].idxmax()
                        obs.loc[ss, "C"] = ss_peaks.loc[i, "C"]
                    elif spec=="hncoca":
                        # Set CAm1 shift to strongest peak in this spin system.         
                        i = ss_peaks.loc[:,"Height"].idxmax()
                        obs.loc[ss, "CAm1"] = ss_peaks.loc[i, "C"]
                    elif spec=="hnca":
                        # Set CA shift to strongest peak in this spin system.         
                        i = ss_peaks.loc[:,"Height"].idxmax()
                        obs.loc[ss, "CA"] = ss_peaks.loc[i, "C"]
                    elif spec=="hncocacb":
                        # Use a simple heuristic to guess if peak is CA or CB:
                        # - If only 1 peak, CA if shift >41 ppm, otherwise CB
                        # - If >1 peak, only keep the two with highest 
                        #   (absolute) intensity. 
                        # - If both >48 ppm, the largest shift is CB. 
                        #   Otherwise, the smallest shift is CB
                        if sum(peaks["SS_name"]==ss)==1:
                            if ss_peaks["C"].item()>41:
                                obs.loc[ss, "CAm1"] = ss_peaks["C"].item()
                            else:
                                obs.loc[ss, "CBm1"] = ss_peaks["C"].item()
                        else:
                            ss_peaks["Abs_height"] = ss_peaks["Height"].abs()
                            # Above line throws a SettingWithCopy warning, 
                            # but I can't seem to fix it
                            
                            ss_peaks = ss_peaks.sort_values(by="Abs_height",
                                                            ascending=False)
                            ss_peaks = ss_peaks.iloc[0:2,:]
                            C_max = ss_peaks["C"].max()
                            C_min = ss_peaks["C"].min()
                            if C_max>48 and C_min>48:
                                obs.loc[ss,"CAm1"] = C_min
                                obs.loc[ss,"CBm1"] = C_max
                            else:
                                obs.loc[ss,"CAm1"] = C_max
                                obs.loc[ss,"CBm1"] = C_min
                    elif spec=="hncacb":
                        # Use a simple heuristic to guess if peak is CA or CB:
                        # - If only 1 peak, CA if shift >41 ppm, otherwise CB
                        # - If >1 peak, only keep the two with highest 
                        #   (absolute) intensity. 
                        # - If strongest peak is 41-48 ppm and >twice height of
                        #   next highest, then it's glycine CA
                        # - Else, if both >48 ppm, the largest shift is CB. 
                        #   Otherwise, the smallest shift is CB
                        if sum(peaks["SS_name"]==ss)==1:
#                            print(ss_peaks)
#                            print(ss_peaks["C"])
#                            print(ss_peaks["C"].item())
                            if ss_peaks["C"].item()>41:
                                obs.loc[ss, "CA"] = ss_peaks["C"].item()
                            else:
                                obs.loc[ss, "CB"] = ss_peaks["C"].item()
                        else:
                            ss_peaks["Abs_height"] = ss_peaks["Height"].abs()
                            # Above line throws a SettingWithCopy warning, 
                            # but I can't seem to fix it
                            ss_peaks = ss_peaks.sort_values(by="Abs_height",
                                                            ascending=False)
                            ss_peaks = ss_peaks.iloc[0:2,:]
                            C_max = ss_peaks["C"].max()
                            C_min = ss_peaks["C"].min()
                            if (ss_peaks.iloc[0,:]["C"]>41 and 
                                ss_peaks.iloc[0,:]["C"]<48 and 
                                ss_peaks.iloc[0,:]["Abs_height"] >
                                2*ss_peaks.iloc[1,:]["Abs_height"]):
                                
                                obs.loc[ss,"CA"] = ss_peaks.iloc[0,:]["C"]
                            elif C_max>48 and C_min>48:
                                obs.loc[ss,"CA"] = C_min
                                obs.loc[ss,"CB"] = C_max
                            else:
                                obs.loc[ss,"CA"] = C_max
                                obs.loc[ss,"CB"] = C_min
                    else:
                        print("Spectrum type %s not recognised" % spec)
                        break
            
        self.obs = obs
        return(self.obs)
            
    def import_obs_shifts(self, filename, filetype, SS_num=False):
        """ Import a chemical shift list
        
        filename: Path to text file containing chemical shifts.
        filetype: Allowed values are "naps", "ccpn", "sparky", "xeasy" or 
            "nmrpipe"
            The "ccpn" option is for importing a Resonance table exported from 
            Analysis v2.x. The "naps" option is for importing an unassigned 
            shift table previously exported from NAPS
        SS_num: If true, will extract the longest number from the SS_name and 
        treat it as the residue number. Without this, it is not possible to get
        the i-1 shifts for each spin system.
            
        """
        # Import from file
        if filetype=="naps":
            obs = pd.read_table(filename)
        elif filetype=="ccpn":
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
#        elif filetype == "nmrstar":
#            tmp = nmrstarlib.read_files(filename)
#            return(tmp)
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
        
    def import_testset_shifts(self, filename, remove_Pro=True, 
                          short_aa_names=True):
        """ Import observed chemical shifts from testset data
        
        This function is intended for use with test data only, and is unlikely 
        to work well on 'real' data.
        """
        #### Import the observed chemical shifts
        obs_long = pd.read_table(filename)
        obs_long = obs_long[["Residue_PDB_seq_code","Residue_label",
                             "Atom_name","Chem_shift_value"]]
        obs_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
        # Convert residue type to single-letter code
        if short_aa_names: 
            obs_long["Res_type"] = obs_long["Res_type"].apply(seq1)
            obs_long["SS_name"] = (obs_long["Res_N"].astype(str) + 
                    obs_long["Res_type"])
            obs_long["SS_name"] = [s.rjust(5) for s in obs_long["SS_name"]]
        else:
            obs_long["SS_name"] = (obs_long["Res_N"].astype(str).rjust(4) + 
                    obs_long["Res_type"])
            obs_long["SS_name"] = [s.rjust(7) for s in obs_long["SS_name"]]
            obs_long["Res_type"] = obs_long["Res_type"].apply(seq1)
        obs_long = obs_long.reindex(columns=["Res_N","Res_type","SS_name",
                                             "Atom_type","Shift"])
        
        # Convert from long to wide
        obs = obs_long.pivot(index="Res_N", columns="Atom_type", 
                             values="Shift")
        
        # Add the other columns back in
        tmp = obs_long[["Res_N","Res_type","SS_name"]]
        tmp = tmp.drop_duplicates(subset="SS_name")
        tmp.index = tmp["Res_N"]
        obs = pd.concat([tmp, obs], axis=1)
        
        # Add HADAMAC information
        hadamac_groups = ["VIA","G","S","T","DN","FHYWC","REKPQML"]
        obs["HADAMAC"]=obs["Res_type"]
        for g in hadamac_groups:
            obs["HADAMAC"] = obs["HADAMAC"].str.replace("["+g+"]", g)
        
        # Make columns for the i-1 observed shifts of C, CA and CB
        obs_m1 = obs[list({"C","CA","CB","HADAMAC"}.intersection(obs.columns))]
        obs_m1.index = obs_m1.index+1
        obs_m1.columns = obs_m1.columns + "m1"
        obs = pd.merge(obs, obs_m1, how="left", left_index=True, 
                       right_index=True)
        
        # Restrict to specific atom types
        atom_set = {"H","N","C","CA","CB","Cm1","CAm1","CBm1","HA","HADAMACm1"}
        obs = obs[["Res_N","Res_type","SS_name"]+
                  list(atom_set.intersection(obs.columns))]
        
        obs.index = obs["SS_name"]
        
        if remove_Pro:
            # Remove prolines, as they wouldn't be observed in a real spectrum
            obs = obs.drop(obs.index[obs["Res_type"].isin(["PRO","P"])]) 
        
        self.obs = obs
        return(self.obs)
    
    def export_obs_shifts(self, filename):
        """ Export a chemical shift list.
        
        The purpose of this function is to make a shift list which the user can 
        manually modify/correct, and then reimport.
        """
        
        df = self.obs.melt(id_vars="SS_name", 
                      value_vars=set(self.obs.columns).intersection({"H","N",
                                    "HA","C","CA","CB","Cm1","CAm1","CBm1"}), 
                      var_name="Atom_type", value_name="Shift")
        df = df.sort_values(by=["SS_name", "Atom_type"])
            
        df.to_csv(filename, sep="\t", float_format="%.3f", index=False)
        
        return(df)
        
        
    
                    

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

#path = "/Users/aph516/GitHub/NAPS/data/P3a_L273R/"
##path = "C:/Users/Alex/GitHub/NAPS/data/P3a_L273R/"
#
#a = NAPS_importer()
#
#roots = a.import_hsqc_peaks(path+"hsqc.txt", "ccpn")
#hncacb = a.import_3d_peaks(path+"hncacb_sparky.txt", "sparky", "hncacb", 
#                         assign_nearest_root=True)
#hnco = a.import_3d_peaks(path+"hnco.txt", "ccpn", "hnco", assign_nearest_root=True)
#hncocacb = a.import_3d_peaks(path+"cbcaconh.txt", "ccpn", "hncocacb", assign_nearest_root=True)
#obs = a.find_shifts_from_peaks()
#df = a.export_obs_shifts(path+"naps_shifts.txt")
#obs2 = a.import_obs_shifts(path+"naps_shifts.txt", "naps")

#%%

#tmp = a.import_obs_shifts(path+"ccpn_resonances.txt", "ccpn", SS_num=True)
#tmp = a.import_obs_shifts(path+"sparky shifts.txt", "sparky", SS_num=True)
#tmp = a.import_obs_shifts(path+"Xeasy shifts.txt", "xeasy", SS_num=True)
#tmp = a.import_obs_shifts(path+"nmrpipe_shifts.tab", "nmrpipe", SS_num=True)
#tmp = a.import_obs_shifts(path+"../testset/A001_bmr4032.str.corr.pdbresno",
#                          "nmrstar", SS_num=False)
#
#a.import_hsqc_ccpn(path+"hsqc HADAMAC.txt")
#a.import_3d_peaks_ccpn(path+"hnco.txt", "hnco")
#a.import_3d_peaks_ccpn(path+"cbcaconh.txt", "hncocacb")
#a.import_3d_peaks_ccpn(path+"hncacb.txt", "hncacb")
##a.guess_peak_atom_types()
#
#
#roots = a.roots
#peaklists = a.peaklists
#
#hncocacb = a.peaklists["hncocacb"]
#
## Make a list of all possible peak assignments for the HNCOCACB
#atoms = ["CAm1","CBm1"]
#assignments = pd.DataFrame(columns=["SS_name","As_ID"]+atoms)
#for SS in a.roots["SS_name"]:
#    peaks = a.peaklists["hncocacb"].loc[a.peaklists["hncocacb"]["SS_name"]==SS, "C"]
#    tmp = find_all_assignments(peaks, atoms)
#    tmp["SS_name"] = SS
#    tmp["As_ID"] = range(tmp.shape[0])
#    assignments = assignments.append(tmp, ignore_index=True)
#
## Check which assignments are plausible
#b = pd.concat([assignments, score_plausibility(assignments)], axis=1)
#
#c = b.loc[b[["Rest","Gly","Ser","Thr"]].apply(max, axis=1)>=0,:]
