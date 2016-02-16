#' #############################################################################
#' 
#' Runs serial tempering for corpus C-6 (a corpus created from the 
#' English Wikipedia), which includes categories 
#'  Leopardus (8)
#'  Lynx (8)
#'  Prionailurus (7)
#'  
#' See also: 
#'  lda_fgs_st_est_M_MCSE_real_data.R
#'  lda_fgs_st_hs_diagnostics.R
#'
#' Versions: 
#'  Jan 23, 2015 - Modified data loading, help(cats)
#'  Nov 08, 2015 - Testing and code cleanup 
#'  Jul 29, 2015 - Initial version 
#' 
#' Example run: 
#'  Rscript lda_fgs_st_hs_C6b.R
#'    
#' #############################################################################

library(ldamcmc)
data("cats") # loads data (for alternatives see load_raw_data.R)
setwd('D:/data/lda-hp-data/') # sets the working directory 


# Initialize global variables  ------------------------------------------

spacing            <- 1
max.iter           <- 54000
max.iter.final     <- 1004000
burn.in            <- 4000
SEED               <- 2015;
tuning.iter        <- 3
save.z             <- 0 # save z samples 
save.beta          <- 0 # save beta samples 
save.theta         <- 1 # save theta samples 
save.st.grid.index <- 0 # save selected h indices 
save.lp            <- 0 # save log posterior
save.hat.ratios    <- 1
save.tilde.ratios  <- 0
verbose            <- 0
K                  <- length(class.labels)
V                  <- length(vocab)
rdata.file         <- paste("fgs-st-hs-", ds.name, "-K", K, "-seed", SEED, "-", 
                            format(Sys.time(), "%Y%b%d%H%M%S"), "-itr", 
                            tuning.iter, ".RData", sep = "")


# Creating h-grids  ----------------------------------------------------

sg.alpha.interval  <- .01
sg.alpha.start     <- .15
sg.alpha.end       <- .35

sg.eta.interval    <- .05
sg.eta.start       <- .8
sg.eta.end         <- 1.

interval           <- .01
alpha.start        <- .065 
alpha.end          <- .465 
eta.start          <- .705
eta.end            <- 1.105 

init.st.grid.index <- 53 # h = (.91, .255)

x.axis             <- seq(sg.alpha.start, sg.alpha.end, by=sg.alpha.interval)
y.axis             <- seq(sg.eta.start, sg.eta.end, by=sg.eta.interval)
st.grid            <- gen_meshgrid(x.axis, y.axis) # generate alpha grid (2-D)

x.axis2            <- seq(alpha.start, alpha.end, by=interval)
y.axis2            <- seq(eta.start, eta.end, by=interval)
h.grid             <- gen_meshgrid(x.axis2, y.axis2) # generate alpha grid (2-D)


# Identifies neighbours for each grid point  -------------------------------

st.grid.nbrs       <- get_grid_neighbors(x.axis, y.axis);


# Serial sampling ----------------------------------------------------------

set.seed(SEED); # sets seed for inference 

ptm                <- proc.time();
num.st.grid        <- ncol(st.grid) # number of labels 
init.st.grid.zetas <- array(1., c(num.st.grid, 1)); 
model              <- lda_fgs_st(K, V, ds$wid + 1, doc.N, h.grid, st.grid, 
                                 st.grid.nbrs, init.st.grid.index, 
                                 init.st.grid.zetas, max.iter, burn.in, spacing, 
                                 tuning.iter, save.z, save.beta, save.theta, 
                                 save.st.grid.index, save.lp, save.hat.ratios, 
                                 save.tilde.ratios, verbose, max.iter.final);
ptm                <- proc.time() - ptm;



# Saves every object into a file -------------------------------------------

save.image(rdata.file)

cat("Execution time (lda_fgs_st) = ", ptm[3], "\n");
