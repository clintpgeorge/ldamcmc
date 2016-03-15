#' #############################################################################
#' 
#' This script computes ratios of the estimates of posterior predictive scores 
#' of the LDA models indexed by default hyperparameters h_{DR}, h_{DA}, and 
#' h_{DG} to the estimate of the posterior predictive score of the empirical 
#' Bayes model, for a given corpus. 
#' 
#' Note: This script is used to generate each row of Table 5 in the paper 
#' "Principled Selection of Hyperparameters in the Latent Dirichlet Allocation 
#' Model".
#' 
#' Versions: 
#'  November 10, 2015 - Initial version  
#' 
#' Example runs: 
#'  Rscript lda_fgs_ppc.R "D:/clintpg/workspace/ldamcmc/data-raw" "felines" "category" .5 .155 1983 3000 1000 1
#'  Rscript lda_fgs_ppc.R "D:/clintpg/workspace/ldamcmc/data-raw" "rec" "category" .46 .09 1983 2000 1000 1
#'  Rscript lda_fgs_ppc.R "D:/clintpg/workspace/ldamcmc/data-raw" "med-christian-baseball" "category" .385 .085 2015 2000 1000 1
#' 
#' #############################################################################

## Loads necessary libraries and sets global variables 

rm(list=ls()); # Removes all objects in the current R state 

library(ldamcmc);

# Handles commandline arguments

args           <- commandArgs(TRUE)
data.dir       <- args[1] 
prefix         <- args[2] 
class.type     <- args[3] 
eb.eta         <- as.numeric(args[4])
eb.alpha       <- as.numeric(args[5]) 
SEED           <- as.numeric(args[6]) 
max.iter       <- as.numeric(args[7]) # the maximum number of Gibbs iterations
burn.in        <- as.numeric(args[8]) 
spacing        <- as.numeric(args[9]) 

cat("\nCommandline Arguments:\n")
cat(paste("\ndata.dir:", data.dir))
cat(paste("\nprefix:", prefix))
cat(paste("\nclass.type:", class.type))
cat(paste("\neb.alpha:", eb.alpha))
cat(paste("\neb.eta:", eb.eta))
cat(paste("\nseed:", SEED))
cat(paste("\nmax.iter:", max.iter))
cat(paste("\nburn.in:", burn.in))
cat(paste("\nspacing:", spacing))
cat("\n")


cat("\nLoading data files...\n")

## Loads data 

setwd(data.dir) # Sets the working dir
vocab.file <- paste(prefix, ".ldac.vocab", sep="")
doc.file <- paste(prefix, ".ldac", sep="")
metadata.file <- paste(prefix, ".csv", sep="")

vocab <- readLines(vocab.file);
V <- length(vocab);
documents <- read_docs(doc.file);
doc.metadata <- read.csv2(metadata.file, header=T, sep=';');

ds <- vectorize_docs(documents)
doc.N <- calc_doc_lengths(documents)
class.labels <- levels(doc.metadata[, class.type])
K <- length(class.labels)



## default hyperparameters 

# Gensim package (2010)

d1.alpha.v <- rep(1/K, K);            
d1.eta <- 1/K;
set.seed(SEED)
d1.lppv <- lda_fgs_lppv(K, V, ds$wid+1, doc.N, d1.alpha.v, d1.eta, max.iter, burn.in, spacing)



# Asuncion et al. (2009) 

d2.alpha.v <- rep(.1, K);
d2.eta <- .1;
set.seed(SEED)
d2.lppv <- lda_fgs_lppv(K, V, ds$wid+1, doc.N, d2.alpha.v, d2.eta, max.iter, burn.in, spacing)



# Griffiths and Steyvers (2004)

d3.alpha.v <- rep(50/K, K);                  
d3.eta <- 0.1;
set.seed(SEED)
d3.lppv <- lda_fgs_lppv(K, V, ds$wid+1, doc.N, d3.alpha.v, d3.eta, max.iter, burn.in, spacing)



## based on the empirical Bayes choice  

eb.alpha.v <- rep(eb.alpha, K);       
set.seed(SEED)
eb.lppv <- lda_fgs_lppv(K, V, ds$wid+1, doc.N, eb.alpha.v, eb.eta, max.iter, burn.in, spacing)



## Saves every object into a file 

rdata.file <- paste(prefix, "-ppc-default-eb-h(", eb.eta, ",", eb.alpha, ")-K", 
                    K, format(Sys.time(), "%Y%b%d%H%M%S"), ".RData", sep="")
save.image(rdata.file)
cat("\nRData is saved to:", rdata.file, "\n")

cat("--------------------------------------------------------------\n", sep="");

cat("EB log(ppv) = ", eb.lppv$lppc, "\n", sep="");
cat("Gensim log(ppv) = ", d1.lppv$lppc, "\n", sep="");
cat("Asuncion et al. (2009) log(ppv) = ", d2.lppv$lppc, "\n", sep="");
cat("Griffiths and Steyvers (2004) log(ppv) = ", d3.lppv$lppc, "\n", sep="");

cat("--------------------------------------------------------------\n", sep="");

cat("\\hat{S}(h_{EB}) = ", exp(eb.lppv$lppc-eb.lppv$lppc), "\n", sep="");
cat("\\hat{S}(h_{DR}) = ", exp(d1.lppv$lppc-eb.lppv$lppc), "\n", sep="");
cat("\\hat{S}(h_{DA}) = ", exp(d2.lppv$lppc-eb.lppv$lppc), "\n", sep="");
cat("\\hat{S}(h_{DG}) = ", exp(d3.lppv$lppc-eb.lppv$lppc), "\n", sep="");


## End 
################################################################################


