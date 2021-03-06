#' #############################################################################
#' 
#' Illustrates LDA Hyperparameter Selection (using serial tempering):   
#' 
#' This script uses a synthetic dataset generated based on the LDA hierarchical 
#' model. This assumes symmetric Dirichlet priors for the LDA model and takes 
#' h = (eta, alpha) = (2, 2).   
#' 
#' 
#' Versions: 
#'  Nov 13, 2015 - Added the script to the package ldamcmc   
#'  May 01, 2015 - Initial version   
#' 
#'   
#' Examples: 
#'  Rscript lda_fgs_st_hs_synth_h22.R "~" 1983 101000 1000 2 1000 20 80 4
#' 
#' #############################################################################

rm(list=ls());
library(ldamcmc); 

## Initialize variables

args           <- commandArgs(TRUE)
data.dir       <- args[1] 
SEED           <- as.numeric(args[2])
max.iter       <- as.numeric(args[3])
burn.in        <- as.numeric(args[4])
K              <- as.numeric(args[5])
D              <- as.numeric(args[6])
V              <- as.numeric(args[7])
doc.size       <- as.numeric(args[8])
tuning.iter    <- as.numeric(args[9])

gen.eta        <- 2
gen.alpha      <- 2

spacing        <- 1
prefix         <- "synth"
save.z         <- 0 # save z samples 
save.beta      <- 0 # save beta samples 
save.theta     <- 0 # save theta samples 
save.h.index   <- 0 # save selected h indices 
save.lp        <- 0 # save log posterior
save.hat.ratios <- 1
save.tilde.est <- 1
verbose        <- 0

setwd(data.dir)
fn.prefix      <- paste("fgs-st-hs-", prefix, "-h(", gen.eta, "-", 
                        gen.alpha, ")-k", K, "-v", V, "-d", D, "-seed", 
                        SEED, "-", format(Sys.time(), "%Y%b%d%H%M%S"), sep="")


## Generates the synthetic corpus 

set.seed(1983);
ds             <- gen_corpus(K, V, D, doc.size, gen.alpha, gen.eta);



# Creating h-grids  ----------------------------------------------------


sg.alpha.interval <- .25
sg.eta.interval   <- .25
sg.alpha.start    <- 1
sg.alpha.end      <- 3
sg.eta.start      <- 1
sg.eta.end        <- 3

interval          <- .15
alpha.start       <- .5
alpha.end         <- 6.5
eta.start         <- .5
eta.end           <- 6.5

gen.st.grid.index <- 41 # h = (2, 2) 

x.axis            <- seq(sg.alpha.start, sg.alpha.end, by=sg.alpha.interval)
y.axis            <- seq(sg.eta.start, sg.eta.end, by=sg.eta.interval)
st.grid           <- gen_meshgrid(x.axis, y.axis) # generate alpha grid (2-D)

x.axis2           <- seq(alpha.start, alpha.end, by=interval)
y.axis2           <- seq(eta.start, eta.end, by=interval)
h.grid            <- gen_meshgrid(x.axis2, y.axis2) # generate alpha grid (2-D)


# Identifies neighbours for each grid point  -------------------------------

st.grid.nbrs   <- get_grid_neighbors(x.axis, y.axis);


# Gibbs sampling ----------------------------------------------------------

ptm <- proc.time();
set.seed(SEED); # sets seed for inference 
num.st.grid <- ncol(st.grid) # number of labels 
init.st.grid.zetas <- array(1., c(num.st.grid, 1)); 
model <- lda_fgs_st(K, V, ds$wid, ds$doc.N, h.grid, st.grid, st.grid.nbrs, 
                    gen.st.grid.index, init.st.grid.zetas, max.iter, burn.in, 
                    spacing, tuning.iter, save.z, save.beta, save.theta, 
                    save.h.index, save.lp, save.hat.ratios, save.tilde.est, 
                    verbose);
ptm <- proc.time() - ptm;
cat("Execution time (lda_fgs_st) = ", ptm[3], "\n");

# Saves every object into a file ------------------------------------------

rdata.file <- paste(fn.prefix, "-itr", tuning.iter, ".RData", sep = "")
save.image(rdata.file);
cat("\nSaves the R Session.\n");

