
# Load my functions (will install all my favourite packages too)
source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R")

# Set up dir structure
if(!dir.exists("data")){dir.create("data")}
if(!dir.exists("data/alignments")){dir.create("data/alignments")}

# Important files
refFname <- "data/test_data/Barley_MorexV3_pseudomolecules_first2e6.fasta"
fai <- fread(paste0(refFname,".fai"),select=1:2,col.names=c("chr","length"))
outLDPRsFname <- "lDPRs.rds"

# Plot ID process as we go?
plotLdprId <- F

# Parameter settings
bw <- 200000 # this is the banding size, ie, how close must alignments be together to count. Narrows/widens the "band" of coinsidered alignments around the diagonal.
ll <- 4000 # lower length limit for alignments
qgap <- 500 # for pointwise interrogation, do it in intervals of ...
dbw <- 10000 # this is the bandwidth for the density approximation based on the length-weighted alignments crossing each interrogation point
windowsize <- 2e6 # Can't self-align a whole chromosome without chaining or other filters that would cost us far too much resolution, so instead, do it in chunks of ...
increment <- windowsize-bw # Ensures adequate overlap between windows
cutoff <- 5 # Density cutoff for what is "dense" enough to be considered as a possible l-DPR. Best set based on density plots.

# Lastz settings
#lastzBin <- system("which lastz",intern = T) # Where is the binary
gffFname <- "data/test_data/Barley_MorexV3_pseudomolecules_first2e6.gff3"
lastzBin <- "/apps/easybuild-2022/easybuild/software/MPI/GCC/11.3.0/OpenMPI/4.1.4/LASTZ/1.04.03/bin/lastz" # Where is the binary
lastzArgs <- "--notransition --step=500 --gapped" # Settings for alignment
