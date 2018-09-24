# Script to run NAPS on a single protein.

#### Define protein specific paramaters ####
# Define directories to use
data.dir = "data/P3a L273R"
output.dir = "output/P3a L273R"
plot.dir = "plots/P3a L273R"

# Control and logging parameters
log.file = paste0(output.dir, "/log.txt")
import_peaklists = TRUE     # For if you have peaklists from different spectra to import
import_shifts.obs = FALSE   # If you already have a saved shifts.obs dataframe (eg. if you've manually edited the dataframe created from peak lists)
match_shifts = TRUE
check_consistency = TRUE
plot_results = TRUE

# Parameters for peak list import
peak.file = list(hsqc=paste0(data.dir, "/hsqc HADAMAC.txt"),
                 hnco=paste0(data.dir, "/hnco.txt"),
                 hncaco=NA,
                 cbcaconh=paste0(data.dir, "/cbcaconh.txt"),
                 hncacb=paste0(data.dir, "/hncacb.txt"))

HNC = list(hnco=c(2,3,1),
           hncaco=NA,
           cbcaconh=c(2,3,1),
           hncacb=c(2,3,1))

hncacb.CA.pos = FALSE

observed.shifts.file = paste0(output.dir, "/shifts_obs.txt")

# Parameters for matching peaks to shiftx2 predictions
shiftx2.file = paste0(data.dir, "/shiftx2.cs")
sparta.file = paste0(data.dir, "/P3a_L273R_sparta.tab")
offset = 208
sequence = "TDYIQLLSEIAKEQGFNITYLDIDELSANGQYQCRAELSTSPITVCHGSGISCGNAQSDAAHNALQYLKIIAERK"
seq.start = 239
use.HADAMAC = TRUE
hadamac.error.rate = 0.01
atom.sd = list(H=0.1711, N=1.1169, HA=0.1231, C=0.5330, CA=0.4412, CB=0.5163, C.m1=0.5330, CA.m1=0.4412, CB.m1=0.5163)

assignment.file = paste0(output.dir, "/assignment.txt")
summary.file = paste0(output.dir, "/summary.txt")

# Parameters for plotting
plot.prefix = ""

#### Run NAPS ####
# Load functions
source("R/NAPS_functions.R")

logf = file(log.file, 'w')

source("NAPS_main.R")

close(logf)