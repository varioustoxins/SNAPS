# Script to run NAPS in batch mode.
# Really only intended for testing purposes, on the 61 protein test set used to evaluate shiftx2

library(tidyr)

#### Functions ####
source("R/NAPS_functions.R")

# read_shiftx2_predictions = function(prediction.file, offset=0, ...) {
#   # Imports a chemical shift prediction file from ShiftX2, in csv format
#   # Offset is added to the residue number, in case it is incorrect in the structure
#   shifts.pred.raw = read.csv(prediction.file, as.is=TRUE)
#   shifts.pred.raw$NUM = shifts.pred.raw$NUM + offset # Correct residue numbering
#   # Rename columns, taking account that sometimes there is only sometimes a chain column
#   if (length(names(shifts.pred.raw))==4)   names(shifts.pred.raw) = c("Res.N","Res.type","Atom.type","Shift")
#   else names(shifts.pred.raw) = c("Chain", "Res.N","Res.type","Atom.type","Shift")
#   shifts.pred = spread(shifts.pred.raw, Atom.type, Shift) # Convert from long to wide
#   
#   # Rejig table to get C-1, CA-1 and CB-1 shifts for each residue
#   tmp = shifts.pred[,names(shifts.pred)%in%c("Res.N","Res.type","C","CA","CB")]
#   tmp$Res.N = tmp$Res.N+1
#   shifts.pred = merge(shifts.pred, tmp, by="Res.N", all.x=TRUE, suffixes=c("",".m1"))
#   
#   shifts.pred$Res.name = paste0(shifts.pred$Res.N, shifts.pred$Res.type)
#   return(shifts.pred)
# }

read_observed_shifts = function(filename) {
  shifts.long = read.delim(filename, sep="\t", stringsAsFactors = FALSE)
  shifts.long = shifts.long[,c(2,4,5,7)]
  names(shifts.long) = c("Res.N","Res.type","Atom.type","Shift")
  
  # Convert residue names to single-letter codes
  shifts.long$Res.type = paste0(substr(shifts.long$Res.type,1,1), tolower(substr(shifts.long$Res.type,2,3)))
  library(seqinr)
  shifts.long$Res.type = a(shifts.long$Res.type)  # Convert residue names to single-letter codes
  
  shifts.obs = spread(shifts.long, Atom.type, Shift)
  
  # Rejig table to get C-1, CA-1 and CB-1 shifts for each residue
  tmp = shifts.obs[,names(shifts.obs)%in%c("Res.N","Res.type","C","CA","CB")]
  tmp$Res.N = tmp$Res.N+1
  shifts.obs = merge(shifts.obs, tmp, by="Res.N", all.x=TRUE, suffixes=c("",".m1"))
  
  shifts.obs$Spin.system = paste0(shifts.obs$Res.N, shifts.obs$Res.type)
  return(shifts.obs)
}

read_sequence = function(seq.file, shift.file) {
  # Returns a sequence contained in seqfile.
  # To get the sequence numbering right, also need to read shiftfile
  
  # Read sequence in as a string
  seq = readLines(seq.file)
  
  # Work out the numbering offset, if any
  shifts.long = read.delim(shift.file, sep="\t", stringsAsFactors = FALSE)
  offset = shifts.long[1,2] - shifts.long[1,3]
  
  seq.df = data.frame(Res.N=(1:nchar(seq))+offset, Res.type=strsplit(seq,"")[[1]])
  seq.df$Res.name = paste0(seq.df$Res.N, seq.df$Res.type)
  return(seq.df)
}
#read_sequence("data/testset/BMRB_seqs/6114.txt","data/testset/simplified_BMRB/6114.txt")

#### Main script ####

# Set parameters that stay constant across datasets
# Define directories to use
data.dir = "data/testset"
output.dir = "output/testset"
plot.dir = "plots/testset"

import_peaklists = FALSE
import_shifts.obs = FALSE
match_shifts = TRUE
check_consistency = TRUE
plot_results = TRUE

offset = 0
use.HADAMAC = FALSE
hadamac.error.rate = 0.01
atom.sd = list(H=0.1711, N=1.1169, HA=0.1231, C=0.5330, CA=0.4412, CB=0.5163, C.m1=0.5330, CA.m1=0.4412, CB.m1=0.5163)

# Read a table of metadata
metadata = read.delim(paste0(data.dir,"/testset.txt"), sep="\t", stringsAsFactors = FALSE, header=FALSE)
names(metadata) = c("ID","PDB","BMRB","Resolution","N.AA")
metadata$Spin.systems = NA  # create columns to store results
metadata$Predictions = NA
metadata$Correct = NA

for(i in 1:length(metadata$ID)){
#for(i in 1){
  print(metadata$ID[i])
  
  # Set dataset specific parameters
  log.file = paste0(output.dir, "/", metadata$ID[i], "_log.txt")
  
  shiftx2.file = paste0(data.dir, "/shiftx2_results/", metadata$ID[i],"_",metadata$PDB[i],".cs")
  offset = 0
  sequence=NA
  use.HADAMAC = TRUE
  hadamac.error.rate = 0.01
  atom.sd = list(H=0.1711, N=1.1169, HA=0.1231, C=0.5330, CA=0.4412, CB=0.5163, C.m1=0.5330, CA.m1=0.4412, CB.m1=0.5163)
  
  assignment.file = paste0(output.dir, "/", metadata$ID[i], "_assignment.txt")
  summary.file = paste0(output.dir, "/", metadata$ID[i], "_summary.txt")
  
  plot.prefix = paste0(metadata$ID[i], "_")
  
  logf = file(log.file, 'w')

  # Read in the observed shifts
  shifts.obs = read_observed_shifts(paste0(data.dir, "/simplified_BMRB/",metadata$BMRB[i],".txt"))
  shifts.obs = shifts.obs[!shifts.obs$Res.type=="P",] # Get rid of any spin systems from proline residues, as these would not be observed in HN-detected spectra
  cat(paste0("Read ", length(shifts.obs[,1]), " observed shifts from file: simplified_BMRB/", metadata$BMRB[i],".txt\n"), file=logf)
  
  # Create HADAMAC groups for the observed data
  if (use.HADAMAC) {
    shifts.obs$HADAMAC = ""
    HADAMAC.groups = c("G","S","T","DN","AVI","WHFYC","REKPQML")
    for (g in HADAMAC.groups) {
      mask = grepl(paste0("[",g,"]"), shifts.obs$Res.type.m1)
      shifts.obs$HADAMAC[mask] = g
    }
    cat(paste0("Created HADAMAC column for the observed shifts\n"), file=logf)
  }
  
  # restrict the predicted shifts to those present in the observed dataset
  seq = read_sequence(paste0(data.dir, "/BMRB_seqs/",metadata$BMRB[i],".txt"),paste0(data.dir, "/simplified_BMRB/",metadata$BMRB[i],".txt"))
  sequence = paste0(seq$Res.type, collapse="")
  seq.start = seq$Res.N[1]
  limit_preds_to_seq = TRUE
  
  # Run the script
  source("NAPS_main.R")
  
  # Add into a single dataframe
  if (i==1){
    assign.df$ID=metadata$ID[i]
    assign.df.all = assign.df
  }
  else {
    assign.df$ID=metadata$ID[i]
    assign.df.all = bind_rows(assign.df.all, assign.df)
  }
  
  if(isOpen(logf))  close(logf)
}

#### Analysis of accuracy ####
# Cysteines in disulphide bridges are often named B in this dataset. We don't need to know this, so change to C
mask = assign.df.all$Res.type=="B"
mask[is.na(mask)] = FALSE
assign.df.all$Res.type[mask] = "C"
assign.df.all$Res.name[mask] = paste0(assign.df.all$Res.N[mask], "C")

#### Work out which predicted assignments are correct ####
# First work out which spin systems have a matching prediction
assign.df.all$SS.in.pred = NA
for (id in unique(assign.df.all$ID)){
  mask = assign.df.all$ID == id
  assign.df.all$SS.in.pred[mask] =  assign.df.all$Spin.system[mask] %in% assign.df.all$Res.name[mask]
}
# Then do the same with predictions
assign.df.all$Pred.in.ss = NA
for (id in unique(assign.df.all$ID)){
  mask = assign.df.all$ID == id
  assign.df.all$Pred.in.ss[mask] =  assign.df.all$Res.name[mask] %in% assign.df.all$Spin.system[mask]
}

# If Spin.system matches Res.name, the assignment is correct
assign.df.all$Correct = (assign.df.all$Res.name == assign.df.all$Spin.system)
# If assigned to a dummy residue, set Correct to NA
assign.df.all$Correct[assign.df.all$Dummy.res] = NA
# Create column saying whether it's correctly assigned and why
assign.df.all$Status=factor(paste(assign.df.all$SS.in.pred, assign.df.all$Correct), 
                         levels=c("FALSE FALSE","TRUE NA","TRUE FALSE", "FALSE NA","TRUE TRUE"),
                         labels=c("Wrongly assigned","Wrongly unassigned","Misassigned","Correctly unassigned","Correctly assigned"))

# It's not obvious what the best way of evaluating the assignment should be.
# I've chosen to count only observed spin systems, and count as correct only if:
# - spin system is assigned to the correct residues
# - spin system is assigned to a dummy residue, and the correct residue isn't in the set of predictions

# Update metadata table
metadata$Correctly.assigned = NA
metadata$Correctly.unassigned = NA
metadata$Wrongly.assigned = NA
metadata$Misassigned = NA
metadata$Wrongly.unassigned = NA
for (i in 1:length(metadata$ID)) {
  mask = (assign.df.all$ID == metadata$ID[i])
  tmp = assign.df.all[mask & !assign.df.all$Dummy.ss,]
  metadata$Spin.systems[i] = sum(!assign.df.all$Dummy.ss[mask])
  metadata$Predictions[i] = sum(!assign.df.all$Dummy.res[mask])
  metadata$Correctly.assigned[i] = sum(tmp$Correct, na.rm=TRUE)
  metadata$Correctly.unassigned[i] = sum(is.na(tmp$Correct) & !tmp$SS.in.pred, na.rm=TRUE)
  metadata$Wrongly.assigned[i] = sum(!tmp$Correct & !tmp$SS.in.pred, na.rm=TRUE)
  metadata$Misassigned[i] = sum(!tmp$Correct & tmp$SS.in.pred, na.rm=TRUE)
  metadata$Wrongly.unassigned[i] = sum(is.na(tmp$Correct) & tmp$SS.in.pred, na.rm=TRUE)
}

# Create plots showing:
# 1. The total % of correctly predicted spin systems
# 2. The % of correctly predicted spin systems, broken down into more detail
# 3. The % of correctly predicted spin systems, ignoring those assigned to dummy residues

p1 = ggplot(data=assign.df.all[!assign.df.all$Dummy.ss,]) + 
  geom_bar(aes(x=ID, fill=(Status %in% c("Correctly unassigned","Correctly assigned"))), position="fill") + 
  #scale_y_continuous(breaks=0.1*(0:10)) + ggtitle("% of spin systems correctly assigned") +
  options + scale_fill_brewer(palette="Set1", name="Correct assignment?")

p2 = ggplot(data=assign.df.all[!assign.df.all$Dummy.ss,]) + 
  geom_bar(aes(x=ID, fill=Status), position="fill") + 
  scale_y_continuous(breaks=0.1*(0:10)) + ggtitle("% of spin systems correctly assigned (broken down into more detail)") +
  options + scale_fill_brewer(palette="Set1", name="Correct assignment?")

p3 = ggplot(data=assign.df.all[!assign.df.all$Dummy.ss & !assign.df.all$Dummy.res,]) + 
  geom_bar(aes(x=ID, fill=(Status=="Correctly assigned")), position="fill") + 
  #scale_y_continuous(breaks=0.1*(0:10)) + ggtitle("% of spin systems correctly assigned (excluding dummy assignments)") +
  options + scale_fill_brewer(palette="Set1", name="Correct assignment?")

plt = plot_grid(p1, p2, p3, align="v", ncol=1, rel_heights=c(1,1,1), axis="lr")
plot2file(plt, "Accuracy summary.pdf", plot.dir, orientation="portrait")

write.table(metadata, paste0(output.dir, "/metadata.txt"), sep="\t")
write.table(assign.df.all, paste0(output.dir, "/all_assignments.txt"), sep="\t")

library(scales)

ggplot(data=assign.df.all[!assign.df.all$Dummy.ss,]) + 
  geom_bar(aes(x=Status, y=(..count..)/sum(..count..), fill=Status)) + scale_y_continuous(labels=percent) + 
  ggtitle("% of spin systems correctly assigned (broken down into more detail)") +
  options + scale_fill_brewer(palette="Set1", name="Correct assignment?") 

proportion_correct = sum(!assign.df.all$Dummy.ss & assign.df.all$Status=="Correctly assigned") / sum(!assign.df.all$Dummy.ss)

