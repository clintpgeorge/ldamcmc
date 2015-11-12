#' #############################################################################
#' 
#' This script computes the estimate of Posterior Predictive Values based on the 
#' LDA Full Gibbs sampler (FGS) and the Augmented Collapsed Gibbs sampler (ACGS)
#' for the data set wt16. 
#' 
#' Last updated on: November 10, 2015 
#' 
#' #############################################################################

## Loads libraries and sets global variables 

rm(list=ls()); # Removes all objects in the current R state 
SEED <- 1983
set.seed(SEED);
options(digits=2)

library(ldamcmc)

prefix         <- "wt16"
class.type     <- "category"
eb.alpha       <- .02
eb.eta         <- .56
max.iter       <- 2000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 1


cat("\nLoading data files...\n")

## Loads data 

data("wt16.vocab")
data("wt16.docs")
data("wt16.docs.metadata")

V <- length(wt16.vocab);
ds <- vectorize_docs(wt16.docs)
doc.N <- calc_doc_lengths(wt16.docs)
num.docs <- nrow(wt16.docs.metadata)

class.labels <- levels(wt16.docs.metadata[, class.type])
K <- length(class.labels)


## Computes Log Posterior Predictive Value based on FGS and ACGS  

set.seed(SEED)
eb.fgs.lppv <- lda_fgs_lppv_R(K, V, eb.alpha, eb.eta, 
                          ds$did+1, ds$wid+1, doc.N, 
                          max.iter, burn.in, spacing)
set.seed(SEED)
eb.acgs.lppv <- lda_acgs_lppv_R(K, V, eb.alpha, eb.eta, 
                                ds$did+1, ds$wid+1, doc.N, 
                                max.iter, burn.in, spacing)


cat("EB log(ppv) FGS = ", eb.fgs.lppv, "\n", sep="");
cat("EB log(ppv) ACGS = ", eb.acgs.lppv, "\n", sep="");


## Saves every object into a file 

rdata.file <- paste(prefix, "-fgs-acgs-ppc-h(", eb.eta, ",", eb.alpha, ")-K", K, 
                    "-", format(Sys.time(), "%Y%b%d%H%M%S"), ".RData", sep="")
save.image(rdata.file)

cat("\nRData is saved to:", rdata.file)
