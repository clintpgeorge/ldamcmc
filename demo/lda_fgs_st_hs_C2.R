#' #############################################################################
#' 
#' Runs the serial tempering chain for corpus C-2 (a corpus created from the 
#' 20Newsgroups dataset), which includes categories 
#'  rec.autos
#'  rec.motorcycles
#'  rec.sport.baseball
#'  rec.sport.hockey
#'  
#' See also: 
#'  lda_fgs_st_est_M_MCSE_real_data.R
#'  lda_fgs_st_hs_diagnostics.R
#'
#' Versions: 
#'  Nov 08, 2015 - Testing and code cleanup 
#'  Jul 29, 2015 - Initial version 
#' 
#' Example run: 
#'  Rscript lda_fgs_st_hs_C2.R
#'    
#' #############################################################################

rm(list=ls());
library(ldast)

# Initialize global variables  ------------------------------------------

data.dir       <- "D:/workspace/ldamcmc/data-raw/rec-b" # the dataset path 
prefix         <- "rec-b"
class.type     <- "category"
spacing        <- 1
max.iter       <- 54000
burn.in        <- 4000
SEED           <- 2015;
tuning.iter    <- 6
save.z         <- 0 # save z samples 
save.beta      <- 0 # save beta samples 
save.theta     <- 0 # save theta samples 
save.st.grid.index <- 0 # save selected h indices 
save.lp        <- 0 # save log posterior
save.hat.ratios <- 1
save.tilde.ratios <- 1
verbose        <- 0

cat("\nCommandline Arguments:\n")
cat(paste("\ndata.dir:", data.dir))
cat(paste("\nprefix:", prefix))
cat(paste("\nclass.type:", class.type))
cat(paste("\nspacing:", spacing))
cat(paste("\nmax.iter:", max.iter))
cat(paste("\nburn.in:", burn.in))
cat(paste("\nseed:", SEED))
cat("\n")

set.seed(SEED); # setting seed 
setwd(data.dir) # Sets the working dir

# Loading data files ------------------------------------------------------

cat("\nLoading data files...\n")

vocab.file     <- paste(prefix, ".ldac.vocab", sep = "");
doc.file       <- paste(prefix, ".ldac", sep = "");
metadata.file  <- paste(prefix, ".csv", sep = "");

vocab          <- readLines(vocab.file);
documents      <- read_docs(doc.file);
doc.metadata   <- read.csv2(metadata.file, header = T, sep = ';');

ds             <- vectorize_docs(documents)
doc.N          <- calc_doc_lengths(documents)
num.docs       <- nrow(doc.metadata)


class.labels   <- levels(doc.metadata[, class.type])
K              <- length(class.labels)
V              <- length(vocab)


fn.prefix      <- paste("fgs-st-hs-", prefix, "-K", K, "-seed", SEED, "-", 
                        format(Sys.time(), "%Y%b%d%H%M%S"), sep="")


# Creating h-grids  ----------------------------------------------------

sg.alpha.interval <- .005
sg.eta.interval   <- .01
sg.alpha.start    <- .05
sg.alpha.end      <- .11
sg.eta.start      <- .42
sg.eta.end        <- .48
init.st.grid.index <- 58 # h = (.44, .09) 

interval          <- .005
alpha.start       <- .005
alpha.end         <- .2
eta.start         <- .25
eta.end           <- .65

x.axis            <- seq(sg.alpha.start, sg.alpha.end, by=sg.alpha.interval)
y.axis            <- seq(sg.eta.start, sg.eta.end, by=sg.eta.interval)
st.grid           <- gen_meshgrid(x.axis, y.axis) # generate alpha grid (2-D)

x.axis2           <- seq(alpha.start, alpha.end, by=interval)
y.axis2           <- seq(eta.start, eta.end, by=interval)
h.grid            <- gen_meshgrid(x.axis2, y.axis2) # generate alpha grid (2-D)


# Identifies neighbours for each grid point  -------------------------------

st.grid.nbrs   <- get_grid_neighbors(x.axis, y.axis);

# Serial sampling ----------------------------------------------------------

ptm <- proc.time();
set.seed(SEED); # sets seed for inference 
num.st.grid <- ncol(st.grid) # number of labels 
init.st.grid.zetas <- array(1., c(num.st.grid, 1)); 
model <- lda_fgs_st(K, V, ds$wid + 1, doc.N, h.grid, st.grid, st.grid.nbrs, 
                    init.st.grid.index, init.st.grid.zetas, max.iter, burn.in, 
                    spacing, tuning.iter, save.z, save.beta, save.theta, 
                    save.st.grid.index, save.lp, save.hat.ratios, save.tilde.ratios, 
                    verbose);
ptm <- proc.time() - ptm;
cat("Execution time (lda_fgs_st) = ", ptm[3], "\n");

# Saves every object into a file -------------------------------------------

rdata.file <- paste(fn.prefix, "-itr", tuning.iter, ".RData", sep = "")
save.image(rdata.file);
cat("\nSaves the R Session.\n");
