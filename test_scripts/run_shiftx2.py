# Script to run shiftx2 on the test set of protein structures
# Can't use batch mode because they may have different pH and temperatures

import nmrstarlib
import collections
from os import system
import pandas as pd
from pathlib import Path

path = Path("/Users/aph516/GitHub/NAPS/")

# Get table which links BMRBs to PDBs
#tab = ascii.read(path/"Data/testset/testset.txt") 
testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, 
                           names=["ID","PDB","BMRB","Resolution","Length"])
testset_df.index = testset_df["ID"]

# Get pH and temperature from the BMRB files
dir_starfiles = nmrstarlib.read_files(path/"data/testset/CS-corrected-testset-addPDBresno")

starfiles_list = list(dir_starfiles)

pars = {}
for sf in starfiles_list:
  for k in sf.keys():
    if isinstance(sf[k], collections.OrderedDict):
      if sf[k]["Saveframe_category"] == "sample_conditions":
        #print sf["data"]
        
        # Set default values
        pH = 6
        temp = 298
        
        for x in sf[k]["loop_0"][1]:
          if x["Variable_type"] == "pH" or x["Variable_type"] == "pH*":
            pH =  x["Variable_value"]
            #print "  pH:", pH 
          elif x["Variable_type"] == "temperature":
            temp =  x["Variable_value"]
            #print "  temp:", temp
        
        pars[sf["data"]] = (pH, temp)
  
# Run ShiftX2 on each PDB file
for i in testset_df["ID"]:
  print(testset_df.loc[i, "ID"], testset_df.loc[i, "PDB"], testset_df.loc[i, "BMRB"])
  pdbfile = (path/"data/testset/PDB-testset-addHydrogens"/
             (testset_df.loc[i, "ID"]+"_"+testset_df.loc[i, "PDB"]+".pdbH"))
  outfile = (path/"data/testset/shiftx2_results"/
             (testset_df.loc[i, "ID"]+"_"+testset_df.loc[i, "PDB"]+".cs"))
  
  system("python /opt/shiftx2-mac/shiftx2.py -i %s -o %s -p %s -t %s" % 
         (pdbfile, outfile, pars[str(testset_df.loc[i, "BMRB"])][0], 
                                 pars[str(testset_df.loc[i, "BMRB"])][1]))