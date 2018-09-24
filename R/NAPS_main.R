# Main script for assigning a protein with NAPS. 
# This should usually be called from NAPS_single.R or NAPS_batch.R, which will contain the definitions of protein-specific parameters

#### Load packages and set other options ####
# Install any missing libraries
libs.req = c("tidyr","dplyr","ggplot2","corrplot","seqinr","directlabels", "clue","cowplot")
libs.missing = !(libs.req %in% installed.packages())
if(any(libs.missing)) install.packages(libs.req[libs.missing])

library(tidyr)  # For gather() and spread()
library(dplyr)  # For one_of()
library(ggplot2)
library(corrplot)
library(seqinr)  # For conversion between 1- and 3-letter amino acid codes
library(clue) # For the Hungarian algorithm used in the match
library(cowplot) # Required for making plots with multiple aligned axes

# Set parameters used in main program
HADAMAC.groups = c("G","S","T","DN","AVI","WHFYC","REKPQML")
atom.names = c("H","N","C","C.m1","CA","CA.m1","CB","CB.m1","HA")

# Graphical parameters
options = theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5), 
                             panel.grid.major.y = element_line(colour="black"),
                             panel.grid.major.x = element_line(colour="grey"),
                             panel.grid.minor.y = element_line(colour="grey", linetype="solid"),
                             panel.grid.minor.x = element_blank())
colours = scale_colour_brewer(palette="Set1") 
fills = scale_fill_brewer(palette="Set1")

#### MAIN PROGRAM ####

#### Import peak lists ####
if (import_peaklists) {
  shifts.obs = read_hsqc(peak.file[["hsqc"]])
  cat(paste("Imported", length(shifts.obs$Spin.system),"spin systems from HSQC:", peak.file[["hsqc"]]),"\n", file=logf)
  
  if (!is.na(peak.file[["hnco"]])) {
    hnco = read_hnco(peak.file[["hnco"]], HNC[["hnco"]])
    shifts.obs = merge(shifts.obs, hnco, by="Spin.system", all.x=TRUE, sort=TRUE)
    cat(paste("Imported", length(hnco$Spin.system),"spin systems from HNCO:", peak.file[["hnco"]]),"\n", file=logf)
  }
  
  if (!is.na(peak.file[["hncaco"]])) {
    hncaco = read_hncaco(peak.file[["hncaco"]], HNC[["hncaco"]])
    shifts.obs = merge(shifts.obs, hncaco, by="Spin.system", all.x=TRUE, sort=TRUE)
    cat(paste("Imported", length(hncaco$Spin.system),"spin systems from HN(CA)CO:", peak.file[["hncaco"]]),"\n", file=logf)
  }
  
  if (!is.na(peak.file[["cbcaconh"]])) {
    cbcaconh = read_cbcaconh(peak.file[["cbcaconh"]], HNC[["cbcaconh"]])
    shifts.obs = merge(shifts.obs, cbcaconh, by="Spin.system", all.x=TRUE, sort=TRUE)
    cat(paste("Imported", length(cbcaconh$Spin.system),"spin systems from CBCA(CO)NH:", peak.file[["cbcaconh"]]),"\n", file=logf)
  }
  
  if (!is.na(peak.file[["hncacb"]])) {
    hncacb = read_hncacb(peak.file[["hncacb"]], HNC[["hncacb"]], CA.pos = hncacb.CA.pos)
    shifts.obs = merge(shifts.obs, hncacb, by="Spin.system", all.x=TRUE, sort=TRUE)
    cat(paste("Imported", length(hncacb$Spin.system),"spin systems from HNCACB:", peak.file[["hncacb"]]),"\n", file=logf)
  }
  
  # Export the observed shifts as a tab-delimited text file
  write.table(shifts.obs, file=observed.shifts.file, sep="\t", row.names=FALSE)
  cat("Wrote observed chemical shifts for ", length(shifts.obs[,1]), " residues to: ", observed.shifts.file, "\n\n", file=logf)
}

if (import_shifts.obs) shifts.obs = read.delim(observed.shifts.file, stringsAsFactors=FALSE)

#### Match observed shifts to the shiftx2 predictions ####
#print("match")
if (match_shifts) {
  shifts.pred = read_shiftx2_predictions(shiftx2.file, offset, sequence, seq.start) 
  #shifts.pred = read_sparta_predictions(sparta.file, offset, sequence, seq.start) 
  cat(paste0("Read shift predictions for ", length(shifts.pred[,1]), " residues from: ", shiftx2.file, "  Numbering offset was ", offset, ".\n"), file=logf)
  
  if(exists("limit_preds_to_seq")) {  # For testing purposes
    if(limit_preds_to_seq) {
      shifts.pred = shifts.pred[shifts.pred$Res.name %in% seq$Res.name,] # Get rid of predictions for residues that aren't in the sample
      tmp = seq[!(seq$Res.N %in% shifts.pred$Res.N),]
      #shifts.pred = merge(seq, shifts.pred, all.x=TRUE) # Make sure there are predictions for all residues in sample
    }
  }
  
  # Limit observed data to backbone atoms
  shifts.obs = shifts.obs[,intersect(c("Spin.system", "Res.N", "Res.type", "Res.type.m1", "HADAMAC", atom.names), names(shifts.obs))]
  
  # Discard shift predictions for which there is no observed data
  atoms.obs = intersect(names(shifts.obs), atom.names)
  
  shifts.pred = shifts.pred[,c("Res.N","Res.name","Res.type","Res.type.m1", atoms.obs)]
  
  scores = calc_score_matrix2(shifts.obs, shifts.pred)
  
  # Penalise residues where the HADAMAC results don't match the predicted residue type
  if (use.HADAMAC) {
    scores = score_HADAMAC(scores, shifts.obs, shifts.pred, mismatch.multiplier=hadamac.error.rate)
    cat("Applied penalty to HADAMAC mismatches.\n", file=logf)
  }
  
  # Create dummy residues to account for prolines, and make the scores matrix square
  tmp = add_dummies(shifts.obs, shifts.pred, scores)
  shifts.obs = tmp[[1]]
  shifts.pred = tmp[[2]]
  scores = tmp[[3]]
  
  # Propose assignment by taking most likely predicted residue for each observed spin system
  assign.df = construct_assignment2(shifts.obs, shifts.pred, scores)
  cat("Generated most likely assignment.\n", file=logf)

  # For each pairing in the assignment, generate 2nd and 3rd choice assignments
  # ?
  
  # Check whether the i, i-1 shifts of sequential residues are consistent
  if(check_consistency) {
  if(any(c("C","CA","CB") %in% names(assign.df) & c("C.m1","CA.m1","CB.m1") %in% names(assign.df))) {
    assign.df = test_consistency(assign.df)
    cat("Checked whether sequential links were consistent.\n", file=logf)
  }
  }
  
  # Export the proposed assignments as a tab-delimited text file
  write.table(assign.df, file=assignment.file, sep="\t", row.names=FALSE)
  cat("Wrote proposed assignment to file: ", assignment.file, "\n", file=logf)
  
  # Make a results table
  summary.df = data.frame(Parameter=c("Spin systems", "Predictions", "Max log score"), 
                          Value=c(sum(!assign.df$Dummy.ss), 
                                  sum(!assign.df$Dummy.res),
                                  max(assign.df$Log.score)))
  write.table(summary.df, file=summary.file, sep="\t", row.names=FALSE)
  cat("Wrote assignment summary to file: ", summary.file, "\n", file=logf)
}

#### Plot some graphs ####
#print("plot")
if (plot_results) {
  # Import the assignment file
  assign.df = read.delim(assignment.file, stringsAsFactors=FALSE)
  
  if("Max.inconsistency.p1" %in% names(assign.df)){
    plot2file(plot_strips3(assign.df), paste0(plot.prefix, "Assignment strip plot.pdf"), plot.dir, size=c(max(11.69,11.69/80*length(assign.df[,1])), 8.27))
  }
  
  #plot2file(plot_score_matrix(scores, assign.df), paste0(plot.prefix, "Normalised score matrix.pdf"), plot.dir)
  
  plot2file(plot_log_score(assign.df), paste0(plot.prefix, "Log probability.pdf"), plot.dir)
  
  plot2file(plot_shift_diff(assign.df, shifts.pred), paste0(plot.prefix, "Shift differences.pdf"), plot.dir)
  
}