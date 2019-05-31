# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 18:58:03 2019

@author: kheyam
"""

import numpy as np
import pandas as pd
from NAPS_importer import NAPS_importer
from NAPS_assigner import NAPS_assigner
from NAPS_analyse import import_testset_metadata, collect_assignment_results, summarise_results
from pathlib import Path
from plotnine import *
import numpy.random as rand

path = Path("..")

testset_df = import_testset_metadata(path)

id_all_carbons = [
        'A002', 'A003', 'A004', 'A005', 'A006', 'A008', 'A009', 'A010', 'A011',
        'A012', 'A013', 'A014', 'A015', 'A016', 'A017', 'A018', 'A019', 'A020',
        'A021', 'A023', 'A025', 'A026', 'A027', 'A028', 'A029', 'A033', 'A035',
        'A036', 'A037', 'A039', 'A043', 'A044', 'A045', 'A049', 'A050', 'A051',
        'A053', 'A059', 'A061', 'A062', 'A066', 'A067', 'A069']

#%%
def remove_random_obs(a, fraction, remove="SS", seed=None):
    """Remove a percentage of observed data for testing purposes.
    
    fraction: the proportion of observations to delete
    remove: can be either "SS" or "carbons"
    """

    if seed is not None:
        rand.seed(seed)
    
    to_keep = (1-fraction)*len(a.obs.index)
    if remove=="SS":
        to_keep = int((1-fraction)*len(a.obs.index))
        selection = list(rand.choice(a.obs.index,to_keep, replace=False))
        a.obs = a.obs.loc[selection,:]
    elif remove=="carbons":
        for atom in {"C","CA","CB","C_m1","CA_m1","CB_m1"}.intersection(a.obs.columns):
            tmp = a.obs.index[~a.obs[atom].isna()]
            to_remove = int(fraction*len(tmp))
            selection = list(rand.choice(tmp,to_remove, replace=False))
            a.obs.loc[selection,atom] = np.NaN

def test_removed_obs(output_dir, id_list, fraction_remove, remove, seed=None, atom_set=None):
    # Create output directory if it doesn't already exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    for id in id_list:
        print(id)
        # Import observed and (simulated) predicted shifts
        importer = NAPS_importer()
        importer.import_testset_shifts(testset_df.loc[id, "obs_file"])
        
        a = NAPS_assigner()
        a.read_config_file(path/"config/config.txt")   
        
        if atom_set is not None:
            a.pars["atom_set"] = atom_set
        
        a.obs = importer.obs
        remove_random_obs(a, fraction_remove, remove, seed)
        
        a.import_pred_shifts(testset_df.loc[id, "preds_file"], "shiftx2")
        
        a.add_dummy_rows()
        a.calc_log_prob_matrix2(sf=1, verbose=False)
        matching = a.find_best_assignments()
        a.make_assign_df(matching, set_assign_df=True)
        a.check_assignment_consistency(threshold=0.1)
        a.assign_df.to_csv(Path(output_dir)/(testset_df.loc[id, "out_name"]+".txt"), 
                           sep="\t", float_format="%.3f", index=False)
        
    assigns = collect_assignment_results(output_dir, testset_df, id_list)
    summary = summarise_results(assigns)
    
    return(assigns, summary)
    
#%% Try deleting spin systems
assigns_dict = {}
summary_dict = {}
rand.seed(0)

for x in [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]:
    print("#### "+str(x))
    assigns_dict[x], summary_dict[x] = test_removed_obs(path/("output/removed_SS/remove_"+str(x)), id_all_carbons, x, "SS")

keys = list(summary_dict.keys())
tmp = pd.DataFrame({"Fraction_removed":keys, 
                    "Pc_correct":[summary_dict[k].loc[0, "Pc_correct"] for k in keys],
                    "High_pc":[summary_dict[k].loc[0, "High_pc"] for k in keys],
                    "Medium_pc":[summary_dict[k].loc[0, "Medium_pc"] for k in keys],
                    "Low_pc":[summary_dict[k].loc[0, "Low_pc"] for k in keys],
                    "Likely wrong_pc":[summary_dict[k].loc[0, "Likely wrong_pc"] for k in keys]})

#### Make a poster figure
plt = ggplot(data=tmp) + geom_point(aes(x="100*Fraction_removed", y="Pc_correct", colour="'a'"))
plt += geom_point(aes(x="100*Fraction_removed", y="High_pc", colour="'b'"))
plt += geom_point(aes(x="100*Fraction_removed", y="Medium_pc", colour="'c'"))
plt += geom_point(aes(x="100*Fraction_removed", y="Low_pc", colour="'d'"))
plt += scale_colour_brewer("qual", palette=6, name="",
                           labels=["Overall accuracy","% High confidence","% Medium confidence","% Low confidence"])
plt = plt + ylab("Percentage") + xlab("% spin systems removed")
plt = plt + scale_y_continuous(breaks=np.arange(0,101,10), limits=(0,100)) 
plt += theme_bw()
plt.save(path/"plots/Poster remove SS accuracy.pdf", height=100, width=100, units="mm")

    
#%% Try deleting carbon shifts
assigns_dict2 = {}
summary_dict2 = {}
rand.seed(1)

for x in [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]:
    print("#### "+str(x))
    assigns_dict2[x], summary_dict2[x] = test_removed_obs(path/("output/removed_carbons/remove_"+str(x)), id_all_carbons, x, "carbons")

keys2 = list(summary_dict2.keys())
tmp2 = pd.DataFrame({"Fraction_removed":keys2, 
                    "Pc_correct":[summary_dict2[k].loc[0, "Pc_correct"] for k in keys2],
                    "High_pc":[summary_dict2[k].loc[0, "High_pc"] for k in keys2],
                    "Medium_pc":[summary_dict2[k].loc[0, "Medium_pc"] for k in keys2],
                    "Low_pc":[summary_dict2[k].loc[0, "Low_pc"] for k in keys2],
                    "Likely wrong_pc":[summary_dict2[k].loc[0, "Likely wrong_pc"] for k in keys2]})

plt = ggplot(data=tmp2) + geom_point(aes(x="100*Fraction_removed", y="Pc_correct", colour="'a'"))
plt += geom_point(aes(x="100*Fraction_removed", y="High_pc", colour="'b'"))
plt += geom_point(aes(x="100*Fraction_removed", y="Medium_pc", colour="'c'"))
plt += geom_point(aes(x="100*Fraction_removed", y="Low_pc", colour="'d'"))
plt += scale_colour_brewer("qual", palette=6, name="",
                           labels=["Overall accuracy","% High confidence","% Medium confidence","% Low confidence"])
plt = plt + ylab("Percentage") + xlab("% carbons removed")
plt = plt + scale_y_continuous(breaks=np.arange(0,101,10), limits=(0,100)) 
plt += theme_bw()
plt.save(path/"plots/Poster remove carbons accuracy.pdf", height=100, width=100, units="mm")
