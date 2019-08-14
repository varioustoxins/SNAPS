# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 21:53:43 2019

@author: Alex
"""
import pandas as pd
from NAPS_importer import NAPS_importer
import argparse
#from pathlib import Path
import logging

parser = argparse.ArgumentParser(description="NAPS (NMR Assignments from Predicted Shifts)")
# Arguments:
# - File containing metadata on peaklists to import
# - Where to write the output
parser.add_argument("input_file", 
                    help="A table containing details of peaklist files.")
parser.add_argument("output_file")
parser.add_argument("-l", "--log_file", default=None)

if True:
    args = parser.parse_args()
else:
    args = parser.parse_args(("../data/peaklists.txt",
                              "../output/test.txt",
                              "-l", "../output/test.log"))
    
# Set up logging
if isinstance(args.log_file, str):
    print(args.log_file)
    logging.basicConfig(filename=args.log_file, filemode='w', level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s")
else:
    logging.basicConfig(level=logging.ERROR)

#%%
importer = NAPS_importer()
    
peaklist_info = pd.read_table(args.input_file, sep="\s+", comment="#")

# Import the root hsqc peaks
root = peaklist_info.loc[peaklist_info["SS_method"]=="root", :]
importer.import_hsqc_peaks(root["filename"][0], root["filetype"][0])
logging.info("Read in %d peaks from %s.", 
             len(importer.roots["SS_name"]), root["filename"][0])

# Import the 3D spectra
for i in peaklist_info.index:
    if peaklist_info.loc[i, "SS_method"]=="root":
        pass
    elif peaklist_info.loc[i, "SS_method"]=="closest_root":
        tmp = importer.import_3d_peaks(peaklist_info.loc[i, "filename"],
                                       peaklist_info.loc[i, "filetype"],
                                       peaklist_info.loc[i, "spectrum"],
                                       assign_nearest_root=True)
        logging.info("Read in %d peaks from %s.", 
             len(tmp["SS_name"]), peaklist_info.loc[i, "filename"])
    elif peaklist_info.loc[i, "SS_method"]=="from_peaklist":
        tmp = importer.import_3d_peaks(peaklist_info.loc[i, "filename"],
                                       peaklist_info.loc[i, "filetype"],
                                       peaklist_info.loc[i, "spectrum"],
                                       assign_nearest_root=False)
        logging.info("Read in %d peaks from %s.", 
             len(tmp["SS_name"]), peaklist_info.loc[i, "filename"])
importer.find_shifts_from_peaks()
logging.info("Generated shift list with %d spin systems.", 
             len(importer.obs["SS_name"]))
importer.export_obs_shifts(args.output_file)
logging.info("Exported chemical shift file in NAPS format to %s.", 
             args.output_file)