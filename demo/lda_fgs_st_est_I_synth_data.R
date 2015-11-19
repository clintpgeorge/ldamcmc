#' #############################################################################
#' 
#' Estimation of a Family of Posterior Expectaions using Serial Tempering: 
#' 
#' This script helps to estimate the posterior probability that two 
#' documents have similar topics, e.g., \code{g(\theta) = I(\l \theta_i - 
#' \theta_j \l < epsilon)}. It supports or tuning for the label parameters as 
#' in Geyer (2011), Section 11.2.4. 
#' 
#' References: 
#'  1. Principled Selection of Hyperparameters in the Latent Dirichlet  
#'     Allocation Model (George and Doss 2015)
#'  2. Importance Sampling, Simulated Tempering, and Umbrella Sampling (Geyer 
#'     2011) 
#' 
#' Versions:
#'  May 30, 2015 - Input variable name changes 
#'  May 28, 2015 - Adding argument handling  
#'  May 23, 2015 - Testing using a new dataset 
#'  April 30, 2015 - Runs the script using Doss's tuning scheme 
#'  March 03, 2015 - Runs the script using the new package "ldast" 
#'  January 19, 2015 - Implemented pseudo-priors tuning in the C++ code 
#'  January 03, 2015 - Rechecking after the x'mas break 
#'  November 04, 2014 - Initial version 
#' 
#' Example runs: 
#'  Rscript lda_fgs_st_est_I_synth_data.R 2015 101000 1000 4
#'  Rscript lda_fgs_st_est_I_synth_data.R 1983 101000 1000 2
#'
#' #############################################################################

rm(list=ls());
library(ldamcmc); 
setwd('H:/lda-hp-data') # the user home directory 

# Initialize variables ----------------------------------------------------

args           <- commandArgs(TRUE)
SEED           <- as.numeric(args[1])
max.iter       <- as.numeric(args[2])
burn.in        <- as.numeric(args[3])
tuning.iter    <- as.numeric(args[4])

gen.alpha      <- .2
gen.eta        <- .4
K              <- 5 # the number of topics
D              <- 20 # the total number of documents to be generated
V              <- 40 # the vocabulary size
doc.size       <- 200
spacing        <- 1
save.z         <- 0 # save z samples 
save.beta      <- 0 # save beta samples 
save.theta     <- 1 # save theta samples 
save.h.index   <- 0 # save selected h indices 
save.lp        <- 0 # save log posterior
save.hat.ratios <- 1
save.tilde.ratios <- 1
verbose        <- 0
fn.prefix      <- paste("fgs-st-exp-synth-k", K, "-v", V, "-d", D, "-ds", 
                        doc.size, "-seed", SEED, "-cfg-02-", 
                        format(Sys.time(), "%Y%b%d%H%M%S"), sep="")

# Generates SYNTHETIC documents -------------------------------------------

set.seed(1983);
ds             <- gen_corpus(K, V, D, doc.size, gen.alpha, gen.eta);


# Creates the h-grid  ----------------------------------------------------

alpha.interval <- .01
alpha.start    <- .1
alpha.end      <- .4
eta.interval   <- .01
eta.start      <- .35
eta.end        <- .45
gen.st.grid.index <- 116 # h = (.4, .2) 

# generates alpha-eta sub-grid (2-D)

x.axis         <- seq(alpha.start, alpha.end, by=alpha.interval)
y.axis         <- seq(eta.start, eta.end, by=eta.interval)
st.grid        <- gen_meshgrid(x.axis, y.axis) 

# generates alpha-eta h-grid (2-D)

x.axis2        <- seq(alpha.start, alpha.end, by=.01)
y.axis2        <- seq(eta.start, eta.end, by=.01)
h.grid         <- gen_meshgrid(x.axis2, y.axis2) 



# Identifies neighbors for each grid point  -------------------------------

st.grid.nbrs <- get_grid_neighbors(x.axis, y.axis);


# Gibbs sampling ----------------------------------------------------------

ptm <- proc.time();
set.seed(SEED); # sets seed for inference 
num.st.grid <- ncol(st.grid) # number of labels 
init.st.grid.zetas <- array(1., c(num.st.grid, 1)); 
model <- lda_fgs_st(K, V, ds$wid, ds$doc.N, h.grid, st.grid, st.grid.nbrs, 
                    gen.st.grid.index, init.st.grid.zetas, max.iter, burn.in, 
                    spacing, tuning.iter, save.z, save.beta, save.theta, 
                    save.h.index, save.lp, save.hat.ratios, save.tilde.ratios, 
                    verbose);
ptm <- proc.time() - ptm;
cat("Execution time (lda_fgs_st) = ", ptm[3], "\n");


# Saves every object into a file ------------------------------------------

rdata.file <- paste(fn.prefix, "-itr", tuning.iter, ".RData", sep = "")
save.image(rdata.file);
cat("\nThe R Session is saved to:", rdata.file, "\n")


