#### Functions for importing peaklists ####
parse_CCPN_assignments = function(assignments, short.aa=TRUE){
  # Take the assignments from a CCPN peak list, and return a data frame of Chain, Res.N, Res.type and Atom.type
  # Note that this only parses fully assigned peaks. Partially assigned, unassigned or ambiguously assigned peaks will return NA
  
  #Create NA vectors of correct length
  Chain = rep(NA, times=length(assignments))
  Res.N = rep(NA, times=length(assignments))
  Res.type = rep(NA, times=length(assignments))
  Atom.type = rep(NA, times=length(assignments))
  
  # Fill in details for unambiguously assigned residues (eg. " A276LeuH")
  assignments = gsub(" ", "", assignments, fixed=TRUE) # Remove spaces from assignment strings
  tmp = grep("[A-Z]*\\d+[A-Za-z]{3}.*",assignments)  # This picks out all rows with the correct format, with or without a chain identifier
  tmp2 = grep("[A-Z]\\d+[A-Za-z]{3}.*",assignments) # This only picks out rows with a chain identifier
  Chain[tmp2] = substr(as.character(assignments[tmp2]),1,1)
  Res.N[tmp] = as.integer(regmatches(assignments[tmp],regexpr("\\d+",assignments[tmp])))
  Res.type[tmp] = regmatches(assignments[tmp],regexpr("\\D\\D\\D",assignments[tmp]))
  Atom.type[tmp] = substr(as.character(assignments[tmp]), regexpr("\\D\\D\\D",assignments[tmp])+3,nchar(as.character(assignments[tmp])))
  
  # Convert amino acid abbreviations to single letter, if desired
  if (short.aa) {
    Res.type[!is.na(Res.type)]=a(Res.type[!is.na(Res.type)])
  }
  
  return(data.frame(Chain, Res.N, Res.type, Atom.type))
}

parse_CCPN_peaks = function(assign){
  # Take assignments from a CCPN peak list, and return a list of spin system identifiers and atom types 
  
  Spin.system = rep("", times=length(assign))
  Atom.type = rep("", times=length(assign))
  
  # Identify peaks with unknown assignments (eg. "{27}H[53]2 or "{93}[197]")
  mask.unassigned = grep("\\{", assign)
  Spin.system[mask.unassigned] = substr(assign[mask.unassigned], regexpr("\\{", assign[mask.unassigned]), regexpr("}", assign[mask.unassigned]))
  Atom.type[mask.unassigned] = substr(assign[mask.unassigned], regexpr("}", assign[mask.unassigned])+1, regexpr("\\[", assign[mask.unassigned])-1)
  
  # Identify peaks with full assignments (eg. "256Ala", "256AlaN", "B256AlaN")
  mask.assigned = grep("^[[:blank:]]*[[:alpha:]]*[[:digit:]]+[[:upper:]][[:lower:]]{2}", assign)
  
  tmp = parse_CCPN_assignments(assign[mask.assigned], short.aa=FALSE)
  
  Spin.system[mask.assigned] = paste0(tmp$Res.N, a(tmp$Res.type))
  Atom.type[mask.assigned] = as.character(tmp$Atom.type)
  
  return(data.frame(Spin.system, Atom.type, stringsAsFactors=FALSE))
}

read_hsqc = function(file) {
  # Take an HSQC peak list (exported from CCPN), and return a list of spin systems, atom types and H,N shifts
  
  hsqc.raw = read.delim(file, stringsAsFactors=FALSE) # Read in the peak list
  
  hsqc = cbind(parse_CCPN_peaks(hsqc.raw$Assign.F1), HADAMAC=NA, H=hsqc.raw$Position.F1, N=hsqc.raw$Position.F2) # Get the spin system names and discard unneeded data
  
  # If the Details column is used, assume it indicates HADAMAC group
  if (!all(hsqc.raw$Details=="None")) {
    hsqc$HADAMAC = hsqc.raw$Details
    hsqc$HADAMAC[!hsqc$HADAMAC %in% HADAMAC.groups] = NA
  }
  
  hsqc = hsqc[hsqc$Atom.type %in% c("H","N",""),-2]   # Restrict to backbone or unassigned atoms, and get rid of Atom.type column
  hsqc = hsqc[order(hsqc$Spin.system),]   # Sort the table
  return(hsqc)
}

read_hnco = function(hnco.file, HNC=c(2,3,1), filter.backbone=TRUE){
  # Convert an HNCO peak list into a list of C-1 shifts linked to spin system identifier
  # The correspondence between H, N, C' and F1, F2, F3 is given by HNC (so if F1 is H, F2 is C' and F3 is N, then HNC=c(1,3,2))
  # If filter.backbone is true, only unassigned peaks and peaks with H or N atom types are kept
  
  hnco.raw = read.delim(hnco.file, stringsAsFactors=FALSE)
  
  # Work out the spin system name
  tmp = parse_CCPN_peaks(hnco.raw[,paste0("Assign.F",HNC[1])])
  hnco.raw$Spin.system = tmp$Spin.system
  
  if (filter.backbone) {
    hnco.raw = hnco.raw[tmp$Atom.type %in% c("","H","N"),] # Keep only backbone or unassigned peaks
  }
  
  hnco = hnco.raw[,c("Spin.system",paste0("Position.F",HNC[3]))] # Just keep columns with spin system name and C' shift
  names(hnco)[2] = "C.m1" # C.m1 == C' of the i-1 residue
  
  return (hnco)
}

read_hncaco = function(hncaco.file, HNC=c(2,3,1), filter.backbone=TRUE){
  # Convert an HN(CA)CO peak list into a list of C shifts linked to spin system identifier
  # The correspondence between H, N, C and F1, F2, F3 is given by HNC (so if F1 is H, F2 is C and F3 is N, then HNC=c(1,3,2))
  # If filter.backbone is true, only unassigned peaks and peaks with H or N atom types are kept
  
  hncaco.raw = read.delim(hncaco.file, stringsAsFactors=FALSE)
  
  # Work out the spin system name
  tmp = parse_CCPN_peaks(hncaco.raw[,paste0("Assign.F",HNC[2])])
  hncaco.raw$Spin.system = tmp$Spin.system
  
  if (filter.backbone) {
    hncaco.raw = hncaco.raw[tmp$Atom.type %in% c("","H","N"),] # Keep only backbone or unassigned peaks
  }
  
  # Keep only the most intesnse peak for each spin system
  hncaco = data.frame(Spin.system = unique(hncaco.raw$Spin.system), C=NA)
  for (i in 1:length(hncaco$Spin.system)) {
    tmp = hncaco.raw[hncaco.raw$Spin.system==hncaco$Spin.system[i],]
    hncaco$C[i] = tmp[which.max(tmp$Height),paste0("Position.F",HNC[3])]
  }
  
  return (hncaco)
}

read_cbcaconh = function(cbcaconh.file, HNC=c(2,3,1), filter.backbone=TRUE){
  # Convert an CBCA(CO)NH peak list into a list of CA and CB shifts linked to spin system identifier
  # The correspondence between H, N, C and F1, F2, F3 is given by HNC (so if F1 is H, F2 is C and F3 is N, then HNC=c(1,3,2))
  # If filter.backbone is true, only unassigned peaks and peaks with H or N atom types are kept
  # Uses a simple heuristic to guess whether peak is CA or CB:
  # - If there's only 1 peak, it's CA if shift is greater than 41 ppm, otherwise it's CB
  # - If there's more than 1 peak, only keep the two with highest (absolute) intensity. 
  # - If both are above 48 ppm, the largest shift is assigned to CB. Otherwise, the smallest shift is CB
  
  cbcaconh.raw = read.delim(cbcaconh.file, stringsAsFactors=FALSE)
  
  # Work out the spin system name
  tmp = parse_CCPN_peaks(cbcaconh.raw[,paste0("Assign.F",HNC[2])])
  cbcaconh.raw$Spin.system = tmp$Spin.system
  
  if (filter.backbone) {
    cbcaconh.raw = cbcaconh.raw[tmp$Atom.type %in% c("","H","N"),] # Keep only backbone or unassigned peaks
  }
  
  # Decide which peaks are CA-1 and CB-1
  cbcaconh = data.frame(Spin.system = unique(cbcaconh.raw$Spin.system), CA.m1=NA, CB.m1=NA)
  for (i in 1:length(cbcaconh$Spin.system)) {
    tmp = cbcaconh.raw[cbcaconh.raw$Spin.system==cbcaconh$Spin.system[i],]
    
    C.col=paste0("Position.F",HNC[3])
    
    if (length(tmp[,1])==1) { # If there's only one peak, assume it's CA unless shift is <40
      if (tmp[1,C.col] > 41) cbcaconh$CA.m1[i] = tmp[1,C.col]
      else cbcaconh$CB.m1[i] = tmp[1,C.col]
      #print(paste(cbcaconh$Spin.system[i],": only one peak"))
    }
    else {
      tmp = tmp[order(abs(tmp$Height), decreasing=TRUE),] # Order by absolute height
      tmp = tmp[1:2,] #Take only the two strongest peaks
      if (all(tmp[,C.col]>48)) { # If both peaks >48 ppm, it's probably serine or threonine.So, smallest shift is CA
        cbcaconh$CA.m1[i] = tmp[which.min(tmp[,C.col]),C.col]
        cbcaconh$CB.m1[i] = tmp[which.max(tmp[,C.col]),C.col]
        #print(paste(cbcaconh$Spin.system[i],": probably S,T"))
      }
      else { # In all other cases, the smallest shift is CB
        cbcaconh$CA.m1[i] = tmp[which.max(tmp[,C.col]),C.col]
        cbcaconh$CB.m1[i] = tmp[which.min(tmp[,C.col]),C.col]
        #print(paste(cbcaconh$Spin.system[i],": other"))
      }
    }
  }
  
  return (cbcaconh)
}

read_hncacb = function(hncacb.file, HNC=c(2,3,1), filter.backbone=TRUE, CA.pos=TRUE){
  # Convert an HNCACB peak list into a list of CA and CB shifts linked to spin system identifier
  # The correspondence between H, N, C and F1, F2, F3 is given by HNC (so if F1 is H, F2 is C and F3 is N, then HNC=c(1,3,2))
  # If filter.backbone is true, only unassigned peaks and peaks with H or N atom types are kept
  # Uses a simple heuristic to guess whether peak is CA or CB:
  # - If there's only 1 peak, it's CA if shift is greater than 41 ppm, otherwise it's CB
  # - If there's more than 1 peak, check if the strongest peak is consistent with glycine (41-48 ppm, and considerably stronger than next nearest peak) 
  # - Otherwise, most intense positive peak is CA and most intense negative peak is CB (adjustable by parameter CA.pos)
  
  hncacb.raw = read.delim(hncacb.file, stringsAsFactors=FALSE)
  
  # Work out the spin system name
  tmp = parse_CCPN_peaks(hncacb.raw[,paste0("Assign.F",HNC[2])])
  hncacb.raw$Spin.system = tmp$Spin.system
  
  if (filter.backbone) {
    hncacb.raw = hncacb.raw[tmp$Atom.type %in% c("","H","N"),] # Keep only backbone or unassigned peaks
  }
  
  # Decide which peaks are CA and CB
  hncacb = data.frame(Spin.system = unique(hncacb.raw$Spin.system), CA=NA, CB=NA)
  for (i in 1:length(hncacb$Spin.system)) {
    tmp = hncacb.raw[hncacb.raw$Spin.system==hncacb$Spin.system[i],]
    
    C.col=paste0("Position.F",HNC[3])
    
    if (length(tmp[,1])==1) { # If there's only one peak, assume it's CA unless shift is <40
      if (tmp[1,C.col] > 41) hncacb$CA[i] = tmp[1,C.col]
      else hncacb$CB[i] = tmp[1,C.col]
      #print(paste(hncacb$Spin.system[i],": only one peak"))
    }
    else {
      tmp = tmp[order(abs(tmp$Height), decreasing=TRUE),] # Order by absolute height
      
      # If strongest peak's shift is consistent with glycine, and has more than twice intensity of next strongest peak, assume glycine.
      if(tmp[1,C.col]<48 & tmp[1,C.col]>41 & abs(tmp$Height[1])>2*abs(tmp$Height[2])) {
        hncacb$CA[i] = tmp[1,C.col]
        #print(paste(hncacb$Spin.system[i],": glycine"))
      }
      else { # In all other cases, the positive peak is CA while the negative is CB (or the  other way round of CA.pos= FALSE)
        h.max = h_max = which.max(tmp$Height)
        h.min = which.min(tmp$Height)
        
        if (CA.pos) {
          hncacb$CA[i] = tmp[h.max,C.col]
          hncacb$CB[i] = tmp[h.min,C.col]
        }
        else {
          hncacb$CA[i] = tmp[h.min,C.col]
          hncacb$CB[i] = tmp[h.max,C.col]
        }
        #print(paste(hncacb$Spin.system[i],": other"))
      }
    }
  }
  
  return (hncacb)
}
#

#### Functions for matching predictions to peaks ####
read_shiftx2_predictions = function(prediction.file, offset=0, sequence=NA, seq.start=1) {
  # Imports a chemical shift prediction file from ShiftX2, in csv format
  # Offset is added to the residue number, in case it is incorrect in the structure
  # The sequence, if provided, is used to fill in any residues for which predictions weren't made
  
  shifts.pred.raw = read.csv(prediction.file, stringsAsFactors=FALSE)
  shifts.pred.raw$NUM = shifts.pred.raw$NUM + offset # Correct residue numbering
  # Rename columns, taking account that sometimes there is only sometimes a chain column
  if (length(names(shifts.pred.raw))==4)   names(shifts.pred.raw) = c("Res.N","Res.type","Atom.type","Shift")
  else names(shifts.pred.raw) = c("Chain", "Res.N","Res.type","Atom.type","Shift")
  shifts.pred = spread(shifts.pred.raw, Atom.type, Shift) # Convert from long to wide
  
  # Use sequence to fill in any missing residues
  # Note: I think this does more harm than good, so have removed it.
  # if(!is.na(sequence)){
  #   tmp = data.frame(Res.N=1:nchar(sequence)+seq.start-1, Res.type=strsplit(sequence, "")[[1]])
  #   shifts.pred = merge(tmp, shifts.pred, by=c("Res.N","Res.type"), all.x=TRUE)
  # }
  
  # Rejig table to get C-1, CA-1 and CB-1 shifts for each residue
  tmp = shifts.pred[,names(shifts.pred)%in%c("Res.N","Res.type","C","CA","CB")]
  tmp$Res.N = tmp$Res.N+1
  shifts.pred = merge(shifts.pred, tmp, by="Res.N", all.x=TRUE, suffixes=c("",".m1"))
  
  shifts.pred$Res.name = paste0(shifts.pred$Res.N, shifts.pred$Res.type)
  return(shifts.pred)
}

read_sparta_predictions = function(prediction.file, offset=0, sequence=NA, seq.start=1) {
  lines = readLines(prediction.file)
  skip.lines = grep("FORMAT", lines) + 1
  
  shifts.pred.raw = read.table(prediction.file, skip=skip.lines, stringsAsFactors=FALSE,
                           col.names=c("Res.N","Res.type","Atom.type","SS_SHIFT","Shift","RC_SHIFT","HM_SHIFT","EF_SHIFT","SIGMA"))
  
  shifts.pred.raw = shifts.pred.raw[,c(1,2,3,5)]
  shifts.pred.raw$Atom.type[shifts.pred.raw$Atom.type=="HN"] = "H"
  #shifts.pred.raw$Res.N = shifts.pred.raw$Res.N + offset
  
  shifts.pred = spread(shifts.pred.raw, Atom.type, Shift)
  shifts.pred$Res.N = shifts.pred$Res.N + offset
  
  # Use sequence to fill in any missing residues
  if(!is.na(sequence)){
    tmp = data.frame(Res.N=1:nchar(sequence)+seq.start-1, Res.type=strsplit(sequence, "")[[1]])
    shifts.pred = merge(tmp, shifts.pred, by=c("Res.N","Res.type"), all.x=TRUE)
  }
  
  # Rejig table to get C-1, CA-1 and CB-1 shifts for each residue
  tmp = shifts.pred[,names(shifts.pred)%in%c("Res.N","Res.type","C","CA","CB")]
  tmp$Res.N = tmp$Res.N+1
  shifts.pred = merge(shifts.pred, tmp, by="Res.N", all.x=TRUE, suffixes=c("",".m1"))
  
  shifts.pred$Res.name = paste0(shifts.pred$Res.N, shifts.pred$Res.type)
  return(shifts.pred)
}
#tmp = read_sparta_predictions("data/P3a L273R/P3a_L273R_sparta.tab", 208)

score_match = function(obs, pred, shifts.pred, atom.sd=list(H=0.1711, N=1.1169, HA=0.1231, 
                                                            C=0.5330, CA=0.4412, CB=0.5163, 
                                                            C.m1=0.5330, CA.m1=0.4412, CB.m1=0.5163)) {
  # Calculate probability of data given prediction. 
  # obs and pred are single rows from the shifts.obs and shifts.pred data frames respectively.
  # Where data is missing, the probability is estimated based on distribution of predicted shifts of this type.
  # k is a scaling constant to adjust the expected likelihood of deviation from the predicted shift.
  
  # Define valid atom types, and the standard deviation between observed and predicted shift (taken from ShiftX2 website)
  atoms =  c("H","N","HA","C","CA","CB","C.m1","CA.m1","CB.m1")
  
  # Convert from wide to long format
  obs = gather(obs, key=Atom.type, value=Shift, one_of(atoms))
  pred = gather(pred, Atom.type, Shift.pred, one_of(atoms))
  
  # Find difference between observed and predicted shift for each atom, and calculate likelihood
  df = merge(pred, obs, by="Atom.type", all.x=TRUE)
  df = df[!is.na(df$Shift.pred),]
  df$diff = df$Shift.pred - df$Shift
  #df$prob = dnorm(df$Shift, mean=df$Shift.pred, sd=as.numeric(atom.sd[df$Atom.type])*k)
  df$prob = 2*pnorm(-abs(df$Shift-df$Shift.pred), sd=as.numeric(atom.sd[df$Atom.type]))
  
  # Deal with cases where observed shift is missing
  for (i in 1:length(df[,1])) {
    if (is.na(df$prob[i])) {
      # Get the list of chemical shifts of the same type from shifts.pred (maybe exclude prolines)
      ref.shifts = shifts.pred[shifts.pred$Res.type!="P", df$Atom.type[i]]
      # Calculate the score for each, then take the mean.
      ref.probs = 2*pnorm(-abs(df$Shift-df$Shift.pred), sd=as.numeric(atom.sd[df$Atom.type]))
      df$prob[i] = mean(ref.probs, na.rm=TRUE)
    }
  }
  
  # If spin system has an observed H shift, it can't be a proline
  if (pred$Res.type[1]=="P" & "H" %in% obs$Atom.type) return(0)
  # Otherwise, calculate product of the individual probabilities
  else return(prod(df$prob))
}

score_match2 = function(shifts.obs, pred, atom.sd=list(H=0.1711, N=1.1169, HA=0.1231, 
                                                       C=0.5330, CA=0.4412, CB=0.5163, 
                                                       C.m1=0.5330, CA.m1=0.4412, CB.m1=0.5163)) {
  # Calculate probability of data given prediction. Calculates a whole set of predictions at once
  # obs is a single row from the shifts.obs data frame.
  # Where data is missing, the probability is estimated based on distribution of predicted shifts of this type.
  
  atoms = intersect(atom.names, intersect(names(shifts.obs), names(pred)))
  
  for (a in atoms) {
    shifts.obs[,paste0("d",a)] = shifts.obs[,a] - pred[,a]
    shifts.obs[,paste0("p",a)] = 2*pnorm(-abs(shifts.obs[,paste0("d",a)]), sd=atom.sd[[a]])
  }
  
  # Deal with NA values
  tmp = shifts.obs[,paste0("p",atoms)]
  tmp[is.na(tmp)] = 0.01
  shifts.obs[,paste0("p",atoms)] = tmp
  
  shifts.obs$Prob = apply(shifts.obs[,paste0("p",atoms)], 1, prod)
  
  # deal with prolines
  if(pred$Res.type=="P") {
    shifts.obs$Prob[!is.na(shifts.obs$H)] = 0
  }
  
  return(shifts.obs$Prob)
}

#tmp = score_match2(shifts.obs, shifts.pred[10,])

calc_score_matrix = function (shifts.obs, shifts.pred) {
  N.obs = length(shifts.obs[,1])
  N.pred = length(shifts.pred[,1])
  scores = matrix(data=NA, N.obs, N.pred, dimnames=list(shifts.obs$Spin.system, shifts.pred$Res.name))
  cat("*** Calculating match scores ***\n", file=logf)
  print("*** Calculating match scores ***")
  for (i in 1:N.obs) {
    tmp = names(shifts.obs)[names(shifts.obs) %in% atom.names]
    cat(paste(i, ": Spin system", shifts.obs$Spin.system[i]),". Atoms observed: ",paste(tmp),"\n", file=logf)
    print(paste(i, ": Spin system", shifts.obs$Spin.system[i]))
    for (j in 1:N.pred) {
      scores[i,j] = score_match(shifts.obs[i,], shifts.pred[j,], shifts.pred, atom.sd)
    }
  }
  
  return(scores)
}

calc_score_matrix2 = function (shifts.obs, shifts.pred) {
  N.obs = length(shifts.obs[,1])
  N.pred = length(shifts.pred[,1])
  scores = matrix(data=NA, N.obs, N.pred, dimnames=list(shifts.obs$Spin.system, shifts.pred$Res.name))
  cat("*** Calculating match scores ***\n", file=logf)
  #print("*** Calculating match scores ***")
  for (i in 1:N.pred) {
    cat(paste(i, ": Residue", shifts.pred$Res.N[i]),"\n", file=logf)
    #print(paste(i, ": Residue", shifts.pred$Res.N[i]))
    
    scores[,i] = score_match2(shifts.obs, shifts.pred[i,], atom.sd)
    
  }
  
  return(scores)
}

score_HADAMAC = function(scores, shifts.obs, shifts.pred, mismatch.multiplier=0) {
  # Use HADAMAC data to narrow down possible matches
  for (g in HADAMAC.groups) {
    matching.obs = (shifts.obs$HADAMAC == g & !is.na(shifts.obs$HADAMAC))  # find the spin systems which match the current group
    matching.preds = (shifts.pred$Res.type.m1 %in% strsplit(g,"")[[1]]) # Find the residue predictions that match current group
    
    # Penalise any scores where the spin system matches the current group, but the residue prediction doesn't 
    scores[matching.obs, !matching.preds] = scores[matching.obs, !matching.preds]*mismatch.multiplier
  }
  return(scores)
}

normalise_scores = function(scores, by="row") {
  # Takes a scores matrix (Rows are observed spin sytems, columns are residue predictions),
  # and scales it so that all rows (or columns) sum to unity. Also removes NaN's and infinities.
  scores[is.na(scores)] = 0
  
  if (by=="row") {
    row.sums = apply(scores, 1, sum, na.rm=TRUE)
    scores.norm = diag(1/row.sums) %*% scores
    row.names(scores.norm) = row.names(scores)
  }
  else if (by=="col") {
    col.sums = apply(scores, 2, sum, na.rm=TRUE)
    scores.norm = scores %*% diag(1/col.sums)
    colnames(scores.norm) = colnames(scores)
  }
  scores.norm[is.nan(scores.norm)] = 0
  scores.norm[is.infinite(scores.norm)] = 0
  return(scores.norm)
}

add_dummies = function(shifts.obs, shifts.pred, scores) {
  # Function that takes N observed spin systems, M predicted residues and an NxM scores matrix, and returns all of length P, P>=max(M,N)
  
  N = length(shifts.obs[,1])
  M = length(shifts.pred[,2])
  Pro = sum(shifts.pred$Res.type=="P")
  
  shifts.obs$Dummy.ss = FALSE
  shifts.pred$Dummy.res=FALSE
  
  # Residues where all the predicted shifts are NA should be set as dummies
  mask = apply(is.na(shifts.pred[,intersect(atom.names, names(shifts.pred))]), 1, all)
  shifts.pred$Dummy.res[mask] = TRUE
  scores[,mask] = 0
  
  # Add rows to shifts.obs and scores to account for prolines in sequence
  if(Pro>0) {
    ss.prolines =  shifts.obs[N+1:Pro,]
    ss.prolines$Spin.system = paste0("proline_",1:Pro)
    ss.prolines$Dummy.ss = TRUE
    shifts.obs = rbind(shifts.obs, ss.prolines)
    
    dummy.prolines = matrix(0, nrow=Pro, ncol=dim(scores)[2])  # Initialise match score as zero.
    dummy.prolines[,shifts.pred$Res.type=="P"] = 1  # Match score is 1 for prolines
    row.names(dummy.prolines) = paste0("Pro_",1:Pro)
    scores = rbind(scores, dummy.prolines)
  }
  
  if (N+Pro > M) {  # If # observations/rows greater than # columns/predictions
    diff = N + Pro - M
    dummy.res = shifts.pred[M+1:diff,]
    dummy.res$Res.name = paste0("dumR_",1:diff)
    dummy.res$Dummy.res = TRUE
    shifts.pred = rbind(shifts.pred, dummy.res)
    
    dummy.cols = matrix(0, nrow=N+Pro, ncol=diff)
    colnames(dummy.cols) = paste0("dumR_",1:diff)
    scores = cbind(scores, dummy.cols)
  }
  else if (N+Pro < M) { # If # observations/rows less than than # columns/predictions
    diff = M - N - Pro
    dummy.ss = shifts.obs[N+Pro+1:diff,]
    dummy.ss$Spin.system = paste0("dumS_",1:diff)
    dummy.ss$Dummy.ss = TRUE
    shifts.obs = rbind(shifts.obs, dummy.ss)
    
    dummy.rows = matrix(0, nrow=diff, ncol=M)
    row.names(dummy.rows) = paste0("dumS_",1:diff)
    scores = rbind(scores, dummy.rows)
  }
  
  return(list(shifts.obs, shifts.pred, scores))
}

#tmp = add_dummies(shifts.obs, shifts.pred, scores)

construct_assignment = function(shifts.obs, shifts.pred, scores) {
  # This takes a spin system list, a prediction list and a scores matrix, and generates a table of proposed assignment
  
  assign.df = shifts.obs
  assign.df$Res.N = shifts.pred$Res.N[apply(scores, 1, which.max)]  # Work out the most likely residue for each spin system
  assign.df$Score = apply(scores, 1, max) # Get the matching score for that residue
  
  #assign.df = assign.df[!(assign.df$CA==0 & assign.df$CB==0 & assign.df$`CA-1`==0 & assign.df$`CB-1`==0),] # Get rid of residues with missing data
  assign.df = merge(shifts.pred[,c("Res.N", "Res.name")], assign.df, by="Res.N", all=TRUE)  # Create rows for residues missing observed data, and get the Res.names
  
  # For residues with multiple predicted spin systems, mark the best one
  assign.df = assign.df[order(assign.df$Res.N, assign.df$Score, decreasing=c(FALSE,TRUE), method="radix"),] # Sort by Residue, then Score
  assign.df$Primary.match = TRUE # Make a column to keep track of whether assignment is the primary match or not
  duplicates = assign.df[duplicated(assign.df$Res.N),]  # Select rows with duplicate residue assignments (these will be the lowest scoring ones)
  duplicates$Res.name = paste0(duplicates$Res.name, "?") # Mark them as less certain
  duplicates$Primary.match = FALSE
  #duplicates$Res.N = NA # This allows duplicates to be sorted to end of list 
  assign.df[duplicated(assign.df$Res.N),] = duplicates # Copy back into main data frame
  
  return(assign.df)
}

construct_assignment2 = function(shifts.obs, shifts.pred, scores) {
  # This takes a spin system list, a prediction list and a scores matrix, and generates a table of proposed assignment
  # Treats it as a matching problem, and uses the Hungarian algorithm to find a maximal probability matching
  
  # Take -log10 of the scores matrix. This sum of the log scores is the log of the total probability of a proposed assignment
  scores.log = -log10(scores)
  
  # Remove infinities and NaNs from the log score matrix
  scores.log[scores.log %in% c(Inf, -Inf, NaN)] = 2*max(scores.log[!(scores.log %in% c(Inf, -Inf, NaN))], na.rm=TRUE)
  
  #Perform the matching
  library(clue)
  matching = solve_LSAP(scores.log)
  
  #Make the assignment dataframe
  assign.df = shifts.obs
  
  assign.df$Res.N = shifts.pred$Res.N[matching]
  assign.df$Res.name = shifts.pred$Res.name[matching]
  assign.df$Res.type = shifts.pred$Res.type[matching]
  
  assign.df$Log.score = scores.log[matrix(c(1:length(assign.df$Res.N),matching), ncol=2)]
  
  # Also get the predicted shifts
  assign.df = merge(assign.df, shifts.pred, by=c("Res.N", "Res.name","Res.type"), all.x=TRUE, suffixes=c("",".pred"))
  
  # Things to do with how good the match is
  # assign.df$Primary.match = TRUE
  # assign.df$Primary.match = !is.na(assign.df$Res.N)
  
  return(assign.df)
}

#tmp = construct_assignment2(shifts.obs, shifts.pred, scores)

generate_slternative_assignments = function(shifts.obs, shifts.pred, scores, assign.df) {
  # For a given assignment, generate the next best assignment for each peak.
  # Do this by repeating the assignment but with the score for particular match set to a very high value.
  
  # Get the matching vector for the best assignment
  
  # For each residue
    # Make a new score matrix with modi

    # redo the assignment
}





score_consistency = function(obs.m1, obs.i, scale.factors=c(1,1,1)) {
  # Function takes two observations, and tests whether their i and i-1 shifts match up
  # Returns the r.m.s deviation between carbon shifts of the i and i-1 residues, each scaled by scale factor
  atoms =  c("H","N","HA","C","CA","CB","C.m1","CA.m1","CB.m1")
  scales = data.frame(Atom.type=c("CA.m1","CB.m1","C.m1"), Scale=scale.factors)
  
  obs.i = gather(obs.i, Atom.type, Shift, one_of(atoms))
  obs.m1 = gather(obs.m1, Atom.type, Shift, one_of(atoms))
  
  obs.i = obs.i[obs.i$Atom.type %in% c("CA.m1","CB.m1","C.m1"),]
  obs.m1 = obs.m1[obs.m1$Atom.type %in% c("CA","CB","C"),]
  obs.m1$Atom.type = paste0(obs.m1$Atom.type, ".m1")
  
  df = merge(obs.i, obs.m1, by="Atom.type", all=TRUE, suffixes=c(".i",".m1"))
  df = merge(df, scales, by="Atom.type", all.x=TRUE)
  df$Diff = df$Shift.i - df$Shift.m1
  df$Diff.scaled = df$Diff / df$Scale
  #print(df)
  #return(1)
  if(all(is.na(df$Diff.scaled))){
    return(NA) 
  }
  else return(sqrt(sum(df$Diff.scaled^2, na.rm=TRUE)))   # Has problem that spin systems with less data will tend to get better scores
}

test_consistency = function(assign.df, threshold=0.5, scale.factor=list(C=1, CA=1, CB=1)) {
  # This function returns assign.df with additional columns giving information on whether the shifts of sequential residues are consistent
  # threshold is the maximum allowed shift difference for sequential residues to count as consistent.
  # scale.factor allows the sensitivity of the threshold to be adjusted for particular atom types.
  
  # Work out which atom types we have sequential information for
  atoms.obs = c("C", "CA", "CB")
  atoms.obs = atoms.obs[atoms.obs %in% names(assign.df) & paste0(atoms.obs,".m1") %in% names(assign.df)]
  
  # Store the dummy rows separately, as they're not needed for the consistency check
  dummy=FALSE
  if(sum(is.na(assign.df$Res.N))>0) {
    dummy=TRUE
    dummy.rows = assign.df[is.na(assign.df$Res.N),]
    assign.df = assign.df[!is.na(assign.df$Res.N),]
  }
  
  tmp = assign.df[,c("Res.N", atoms.obs, paste0(atoms.obs,".m1"))] # Create temporary data frame to calculate shift differences

  # Get the 'i' chemical shifts for the preceeding ('_m1') residue
  tmp2 = assign.df[,c("Res.N", atoms.obs)]
  tmp2$Res.N = tmp2$Res.N + 1
  tmp = merge(tmp, tmp2, by="Res.N", all.x=TRUE, suffixes=c("","_m1"))
  
  # Get the 'i-1' chemical shifts for the succeeding ('_p1') residue
  tmp2 = assign.df[,c("Res.N", paste0(atoms.obs,".m1"))]
  tmp2$Res.N = tmp2$Res.N - 1
  tmp = merge(tmp, tmp2, by="Res.N", all.x=TRUE, suffixes=c("","_p1"))
  
  # Calculate the sequential differences
  for (a in atoms.obs) {
    c1 = paste0(a,".m1")
    tmp[,paste0("d",a,".m1")] = (tmp[,paste0(a,".m1")] - tmp[,paste0(a,"_m1")]) * scale.factor[[a]]
    tmp[,paste0("d",a,".p1")] = (tmp[,a] - tmp[,paste0(a,".m1_p1")]) * scale.factor[[a]]
  }
  #assign.df$Total.inconsistency.m1 = apply(abs(assign.df[,paste0("d", atoms.obs, ".m1")]), 1, sum, na.rm=TRUE)
  #assign.df$Total.inconsistency.p1 = apply(abs(assign.df[,paste0("d", atoms.obs, ".p1")]), 1, sum, na.rm=TRUE)
  tmp$Seq.links.m1 = apply(!is.na(tmp[,paste0("d", atoms.obs, ".m1")]), 1, sum, na.rm=TRUE)
  tmp$Seq.links.p1 = apply(!is.na(tmp[,paste0("d", atoms.obs, ".p1")]), 1, sum, na.rm=TRUE)
  tmp$Max.inconsistency.m1 = apply(abs(tmp[,paste0("d", atoms.obs, ".m1")]), 1, max, na.rm=TRUE)
  tmp$Max.inconsistency.p1 = apply(abs(tmp[,paste0("d", atoms.obs, ".p1")]), 1, max, na.rm=TRUE)
  tmp$Max.inconsistency.m1[tmp$Seq.links.m1==0] = NA
  tmp$Max.inconsistency.p1[tmp$Seq.links.p1==0] = NA
  
  assign.df = merge(assign.df, tmp[,c("Res.N","Max.inconsistency.m1","Max.inconsistency.p1","Seq.links.m1","Seq.links.p1")],
                    by="Res.N", all.x=TRUE)
  
  # Put the dummy rows back
  if(dummy) {
    #dummy.rows[,setdiff(names(assign.df), names(dummy.rows))] = NA # You have to create some extra columns in dummy.rows
    assign.df = bind_rows(assign.df, dummy.rows)
    dummy=FALSE
  }

  # Count the number of consistent and inconsistent links for each residue
  assign.df$Consistent.links = 0
  assign.df$Consistent.links = assign.df$Consistent.links + as.numeric((assign.df$Seq.links.m1>0) & (assign.df$Max.inconsistency.m1<threshold))
  assign.df$Consistent.links = assign.df$Consistent.links + as.numeric((assign.df$Seq.links.p1>0) & (assign.df$Max.inconsistency.p1<threshold))
  assign.df$Inconsistent.links = 0
  assign.df$Inconsistent.links = assign.df$Inconsistent.links + as.numeric((assign.df$Seq.links.m1>0) & (assign.df$Max.inconsistency.m1>threshold))
  assign.df$Inconsistent.links = assign.df$Inconsistent.links + as.numeric((assign.df$Seq.links.p1>0) & (assign.df$Max.inconsistency.p1>threshold))
  
  return(assign.df)
}

#tmp = test_consistency(assign.df)

#### Plotting functions ####
plot_strips = function(assign.df, atoms=NA) {
  # This takes a table of proposed assignments, and plots strips so that sequential links can be seen.
  
  # Limits which atom types are plotted to onlth those with both i and i-1 measurements
  if(is.na(atoms)){
    atoms = c("C", "CA", "CB")
    atoms = atoms[atoms %in% names(assign.df) & paste0(atoms,".m1") %in% names(assign.df)]
  }
  plot.df = gather(assign.df, Atom.type, Shift, one_of(c("H", "N", "HA", "CA", "CB", "C", "CA.m1", "CB.m1", "C.m1"))) # Convert from wide to long (one atom per row)
  # Create column saying if atom is i or i-1
  plot.df$i = 0
  plot.df$i[plot.df$Atom.type %in% c("CA.m1","CB.m1","C.m1")] = -1
  # Correct atom names for i-1 atoms
  plot.df$Atom.type[plot.df$Atom.type=="CA.m1"] = "CA"
  plot.df$Atom.type[plot.df$Atom.type=="CB.m1"] = "CB"
  plot.df$Atom.type[plot.df$Atom.type=="C.m1"] = "C"
  
  plot.df$Group = plot.df$Res.N + plot.df$i # Define a column that will link i peak in strip j to i-1 peak in strip j+1
  plot.df = plot.df[plot.df$Atom.type %in% atoms,] # Limit to atom types we're interested in
  
  # Make the plot
  plt = ggplot(data=plot.df, aes(x=paste0(sprintf("%4s", Res.name), " (", Spin.system, ")"))) + 
    #geom_point(aes(y=Shift, colour=factor(i), shape=endsWith(Res.name, "?")), size=2) + scale_shape_manual(name="Duplicate assignment?", values=c(19,1)) +
    geom_point(aes(y=Shift, colour=factor(i), shape=!(Dummy.res | Dummy.ss)), size=2) + scale_shape_manual(name="Non-dummy assignment?", values=c(1,19)) +
    geom_line(data=plot.df[!plot.df$Dummy.res,], aes(y=Shift, group=Group)) + 
    geom_line(aes(y=Shift, group=Spin.system), linetype="dashed") +
    #geom_text(data=plot.df[plot.df$i==0,], aes(y=Shift, label=Spin.system), vjust=0) +
    xlab("Residue (Spin system)") + ylab("Chemical shift") + ggtitle("Proposed assignments") + 
    facet_grid(Atom.type ~ ., scales="free_y") + scale_y_reverse()
  
  plt = plt + geom_line(data=plot.df[!plot.df$Dummy.res,], aes(y=Shift, group=Group))

  plt = plt + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                 legend.position="top",
                                 panel.grid.major.y = element_blank(),
                                 panel.grid.major.x = element_line(colour="grey"),
                                 panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank()) + 
    scale_colour_brewer(palette="Set1", name="Peak type", labels=c("i-1","i"))
  
  return(plt)
}

plot_strips2 = function(assign.df, atoms=NA) {
  p1 = plot_strips(assign.df, atoms) + theme(legend.position = "bottom", plot.title=element_blank())
  p2 = ggplot(data=assign.df, aes(x=sprintf("%4s",Res.name))) + #scale_x_discrete(breaks=NULL) +
    geom_bar(aes(y=Max.inconsistency.p1), stat="identity") + geom_bar(aes(y=-Max.inconsistency.m1), stat="identity") +
    ylab("Max shift difference") + options + fills + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
  plt = plot_grid(p2, p1, align="v", ncol=1, rel_heights=c(1,3), axis="lr")
  return(plt)
}

plot_strips3 = function(assign.df, atoms=NA) {
  p1 = ggplot(data=assign.df, aes(x=sprintf("%4s",Res.name))) +
    geom_bar(aes(y=Log.score), stat="identity") + ylab("Log score") + 
    ylim(0,max(assign.df$Log.score[!assign.df$Dummy.res])) +
    options + fills + scale_x_discrete(position="top") + theme(axis.title.x = element_blank())
  
  if("Max.inconsistency.p1" %in% names(assign.df)){
    p2 = ggplot(data=assign.df, aes(x=sprintf("%4s",Res.name))) + 
      geom_bar(aes(y=Max.inconsistency.p1), stat="identity") + geom_bar(aes(y=-Max.inconsistency.m1), stat="identity") +
      ylab("Max shift difference") + options + fills + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  p3 = plot_strips(assign.df, atoms) + theme(legend.position = "bottom", plot.title=element_blank())
  
  if("Max.inconsistency.p1" %in% names(assign.df)){
    plt = plot_grid(p1, p2, p3, align="v", ncol=1, rel_heights=c(1,1,4), axis="lr")
  }
  else {
    plt = plot_grid(p1, p3, align="v", ncol=1, rel_heights=c(1,1,4), axis="lr")
  }
  return(plt)
}

plot_score_matrix = function(scores, assign.df) {
  # This makes a correlation plot of the scores matrix, sorted according to a proposed assignment
  scores[is.na(scores)] = 0
  scores[assign.df$Dummy.ss,] = 0
  scores[, assign.df$Dummy.res] = 0
  tmp = match(assign.df$Spin.system, row.names(scores))
  scores = scores[tmp[!is.na(tmp)],]
  scores.norm = normalise_scores(scores, by="row")
  
  corrplot(scores.norm, is.corr=FALSE, title="Normalised score matrix")
  
  return(1)
}

plot_shift_diff = function(assign.df, shifts.pred) {
  # Plots the difference between observed and predicted shifts for each spin system of a proposed assignment
  atoms = c("H", "N", "HA", "CA", "CB", "C", "CA.m1", "CB.m1", "C.m1")
  
  # Convert from wide to long
  assign.df = assign.df[,names(assign.df) %in% c("Res.N", "Res.name", "Atom.type","Res.type", "Spin.system", atoms)]
  assign.df = gather(assign.df, Atom.type, Shift, one_of(atoms))
  shifts.pred = gather(shifts.pred, Atom.type, Shift, one_of(atoms))
  
  df = merge(assign.df, shifts.pred, by=c("Res.N", "Res.name", "Atom.type"), all.x=TRUE, suffixes=c(".obs",".pred"))
  df$Diff = df$Shift.obs - df$Shift.pred
  
  plt = ggplot(data=df) + geom_bar(aes(x=paste0(sprintf("%4s", Res.name), " (", Spin.system, ")"), y=Diff), stat="identity") + 
    xlab("Residue (Spin system)") + ylab("Chemical shift difference") + 
    ggtitle("Difference between observed and predicted chemical shifts") +
    facet_grid(Atom.type ~ ., scales="free_y") + options
  
  return(plt)
}

plot_log_score = function (assign.df) {
  plt = ggplot(data=assign.df) + geom_bar(aes(x=paste0(sprintf("%4s", Res.name), " (", Spin.system, ")"), y=Log.score), stat="identity") +
    ggtitle("Residue-by-residue probability of assignment") + xlab("Residue number (spin system)") + ylab("-log(Probability)") + options
  
  return(plt)
}

plot_inconsistency = function(assign.df) {
  plt = ggplot(data=assign.df, aes(x=paste0(sprintf("%4s", Res.name), " (", Spin.system, ")"))) +
    geom_bar(aes(y=Max.inconsistency.p1), stat="identity") + geom_bar(aes(y=-Max.inconsistency.m1), stat="identity") +
    ggtitle("Maximum shift difference with sequential residues") + xlab("Residue number (spin system)") + ylab("Max shift difference") + options + fills
  
  return(plt)
}

plot_hsqc = function(shifts.obs, labels=FALSE) {
  # Plot the HSQC with spin system labels
  plt = ggplot(data=shifts.obs) + geom_point(aes(x=H, y=N)) + scale_x_reverse() + scale_y_reverse() + 
    theme_bw()
  if (labels) plt = plt + geom_text(aes(x=H, y=N, label=Spin.system), vjust=0, hjust=0)
  return(plt)
}

plot_hsqc_probs = function(shifts.obs, shifts.pred, scores, Res.id) {
  plt = plot_hsqc(shifts.obs)
  scores.red = scores[,colnames(scores)==Res.id]
  scores.red.norm = scores.red/sum(scores.red, na.rm=TRUE)
  tmp = cbind(shifts.obs, Score=scores.red.norm)
  plt = plt + geom_point(data=shifts.pred[shifts.pred$Res.name==Res.id,], aes(x=H, y=N), colour="blue")
  plt = plt + geom_text(data=tmp, aes(x=H, y=N, label=sprintf("%.1f %%", Score*100)), hjust=0, vjust=0)
  return(plt)
}

#plot_hsqc_probs(shifts.obs, shifts.pred, scores, "241Y")

plot_hsqc_contours = function(shifts.obs, shifts.pred, Res.id) {
  H.pred = shifts.pred$H[shifts.pred$Res.name==Res.id]
  N.pred = shifts.pred$N[shifts.pred$Res.name==Res.id]
  plt = plot_hsqc(shifts.obs)
  N=200
  tmp = data.frame(H=seq(min(shifts.obs$H), max(shifts.obs$H), length.out=N),
                   N=rep(seq(min(shifts.obs$N), max(shifts.obs$N), length.out=N), each=N))
  tmp$Prob = sqrt(+(H.pred-tmp$H)^2/0.1711^2 + (N.pred-tmp$N)^2/1.1169^2)
  plt = plt + geom_contour(data=tmp, aes(x=H, y=N, z=Prob), binwidth=1)
  return(plt)
}

plot2file = function (plt, file, plot.dir, size="A4", orientation="landscape") {
  if (size=="A4") {
    h = 8.27
    w = 11.69
  }
  else if (size=="A5") {
    h = 11.69/2
    w = 8.27
  }
  else if (length(size)==2) {
    w = size[1]
    h = size[2]
  }
  
  if (orientation=="portrait") {
    tmp = h
    h = w
    w = tmp
  }
  
  pdf(paste0(plot.dir, "/", file), height=h, width=w, useDingbats=FALSE)
  print(
    plt
  )
  dev.off()
}

#### Output functions ####
write_xeasy = function(assign.df, filename) {
  # Writes an xeasy file suitable for import into CCPN
  # Columns needed: atom.ID, shift, SD, atom.type, Res.N
  # To do: Need to account for some atoms only being present as i-1, and check whether atom names need changing
  
  # Convert from wide to long, and discard unneeded columns
  assign.df2 = gather(assign.df, key="Atom.type", value="Shift", intersect(c("H","N","C","CA","CB","HA"), names(assign.df)))
  xeasy.df = assign.df2[order(assign.df2$Res.N),c("Shift","Atom.type","Res.N")]
  xeasy.df = data.frame(ID=1:length(xeasy.df[,1]), Shift=xeasy.df$Shift, SD=0, Atom.type=xeasy.df$Atom.type, Res.N=xeasy.df$Res.N)
  
  xeasy.df = xeasy.df[!is.na(xeasy.df$Res.N),]  # Exclude dummy residues
  xeasy.df$Shift[is.na(xeasy.df$Shift)] = 999 # Replace NA shifts with 999
  
  write.table(xeasy.df, filename, quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(xeasy.df)
}

write_sparky = function(assign.df, filename) {
  # Writes an sparky file suitable for import into CCPN
  # Columns needed: Res.name, Atom.type, Element, Shift, SD, # assignment
  
  # Convert from wide to long, and discard unneeded columns
  assign.df2 = gather(assign.df, key="Atom.type", value="Shift", intersect(c("H","N","C","CA","CB","HA"), names(assign.df)))
  sparky.df = assign.df2[order(assign.df2$Res.N),c("Shift","Atom.type","Res.N", "Res.type")]
  sparky.df = data.frame(Group = paste0(sparky.df$Res.type, sparky.df$Res.N), Atom=sparky.df$Atom.type, Nuc="", Shift=sparky.df$Shift, Sdev=0, Assignments=1, stringsAsFactors = FALSE)
  sparky.df$Nuc[sparky.df$Atom %in% c("H","HA")] = "1H"
  sparky.df$Nuc[sparky.df$Atom=="N"] = "15N"
  sparky.df$Nuc[sparky.df$Atom %in% c("C","CA","CB")] = "13C"
  
  sparky.df = sparky.df[!is.na(sparky.df$Group),]  # Exclude dummy residues
  sparky.df = sparky.df[!is.na(sparky.df$Shift),]  # Exclude residues with NA shifts
  
  write.table(sparky.df, filename, quote=FALSE, row.names=FALSE, col.names=TRUE)
  return(sparky.df)
}

#tmp = write_sparky(assign.df, paste0(output.dir,"/sparky out.txt"))

