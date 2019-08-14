# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 11:39:38 2018

@author: kheyam
"""

assigns1, summary1 = check_assignment_accuracy(path/"output/testset/", prefix="alt_", ranks=[1], N=10)
assigns2, summary2 = check_assignment_accuracy(path/"output/testset/", prefix="alt_", ranks=[2], N=10)
assigns3, summary3 = check_assignment_accuracy(path/"output/testset/", prefix="alt_", ranks=[3], N=10)

assigns1["Rank"]=1
assigns2["Rank"]=2
assigns3["Rank"]=3

assigns_all = pd.concat([assigns1, assigns2, assigns3])

#%%
tmp = assigns_all.loc[assigns_all["ID"]=="A004",:]
plt = ggplot(aes(x="SS_name", y="Log_prob"), data=tmp)
plt = plt + geom_point(aes(colour="Status"))
#plt = plt + geom_line(aes(group="Status"), data=assigns_all.loc[assigns_all["ID"]=="A001" & assigns_all["Status"]=="Correctly assigned",:])
plt = plt + theme(axis_text_x = element_text(angle=90)) #+ ylim(-50, 0) 
plt

#%%
tmp =assigns_all.loc[assigns_all["ID"]=="A002",:]
tmp["Rank"] = tmp["Rank"].astype(str)
plt = ggplot(data=tmp)
plt = plt + geom_boxplot(aes(x="Rank", y="Log_prob", fill="Status"), position="dodge")+ theme(legend_position="top")

#%%
runfile('/Users/aph516/GitHub/NAPS/python/NAPS.py', wdir='/Users/aph516/GitHub/NAPS/python',
        args="peaks /Users/aph516/GitHub/NAPS/data/peaklists.txt "+
        "/Users/aph516/GitHub/NAPS/output/test.txt "+
        "-l /Users/aph516/GitHub/NAPS/output/log.txt")