#' #############################################################################
#' 
#' Estimate Posterior Predictive Values
#' 
#' This script computes the estimate of Posterior Predictive Values based on the 
#' LDA Full Gibbs sampler (FGS) and the Augmented Collapsed Gibbs sampler (ACGS)
#' for the data set wt16. 
#' 
#' Versions: 
#'  November 10, 2015 - Major changes 
#'  November 23, 2015 - Data loading changes 
#' 
#' Example: 
#'  Rscript lda_fgs_ppc_new.R
#' 
#' #############################################################################

## Loads libraries and data, and sets global variables -------------------------

library(ldamcmc)

data("wt16") # Loads data 

eb.alpha           <- .02
eb.eta             <- .56
max.iter           <- 2000 # the maximum number of Gibbs iterations
burn.in            <- 1000
spacing            <- 1
SEED               <- 1983
K                  <- 10 # length(class.labels)
V                  <- length(vocab);
rdata.file         <- paste(ds.name, "-fgs-ppc-h(", eb.eta, ",", eb.alpha, 
                            ")-K", K, "-", format(Sys.time(), "%Y%b%d%H%M%S"), 
                            ".RData", sep = "")


## Computes Log Posterior Predictive Value based on FGS and ACGS ---------------  

set.seed(SEED)
alpha.v            <- rep(eb.alpha, K)
eb.fgs.lppv        <- lda_fgs_lppv(K, V, ds$wid + 1, doc.N, alpha.v, eb.eta, max.iter, burn.in, spacing)

cat("EB log(ppv) FGS = ", eb.fgs.lppv$lppc, "\n", sep = "");


## Saves every object into a file ---------------------------------------------

save.image(rdata.file)
cat("\nRData is saved to:", rdata.file)
