# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 19:08:55 2018

@author: Alex
"""

import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1

class NAPS_importer:
    # Attributes
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
        hsqc = hsqc[["Assign F1","Position F1","Position F2"]]
        hsqc.columns = ["SS_name","H","N"]      
        hsqc.index = hsqc["SS_name"]
        
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
        if isinstance(self.peaks, pd.DataFrame):
            self.peaks = self.peaks.append(peaks, ignore_index=True)
        else:
            self.peaks = peaks
            
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
                    
                    

#%%

#### Test code

a = NAPS_importer()
a.import_hsqc_ccpn("~/GitHub/NAPS/data/P3a_L273R/hsqc HADAMAC.txt")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/hnco.txt", "hnco")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/cbcaconh.txt", "hncocacb")
a.import_3d_peaks_ccpn("~/GitHub/NAPS/data/P3a_L273R/hncacb.txt", "hncacb")
a.guess_peak_atom_types()

roots = a.roots
peaks = a.peaks
shifts = a.shifts
