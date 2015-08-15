#' #############################################################################
#' This script runs the Augmented Collapsed Gibbs sampler and the full Gibbs 
#' sampler of LDA on synthetic datasets 
#' 
#' 
#' 
#' Last modified on: Oct 30, 2014 
#' #############################################################################

## Loads packages 

library(ldamcmc); 
set.seed(1983)

## Initialize variables

prefix         <- "fgs-acgs"
K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 11000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 40
lambda.hat     <- 80
gen.eta        <- 7
gen.alpha      <- 3
gen.alpha.v    <- array(gen.alpha, c(1, K)); # symmetric Dirichlet
gen.eta.v      <- array(gen.eta, c(1, V));   # symmetric Dirichlet
store.z        <- 1                          # store z samples ? 
store.beta     <- 0                          # store beta samples ? 
store.theta    <- 0                          # store theta samples ? 
store.lp       <- 0                          # store log posterior for each iteration 

file.prefix    <- paste(prefix, "-alpha", gen.alpha, "-eta", gen.eta, sep = "")
rdata.file     <- paste(file.prefix, ".RData", sep = "")[1]



## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

# help(gen_synth_corpus)
ds             <- gen_synth_corpus(D, lambda.hat, gen.alpha.v, beta.m)

## The full Gibbs sampling

# help(lda_fgs)
ptm            <- proc.time();
fg.mdl         <- lda_fgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, 
                          max.iter, burn.in, spacing, store.z, store.beta, 
                          store.theta, store.lp);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

# help(lda_acgs)
ptm            <- proc.time();
cg.mdl         <- lda_acgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, 
                           max.iter, burn.in, spacing, store.z, store.beta, 
                           store.theta, store.lp);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Saves all objects into a file

save.image(rdata.file)



