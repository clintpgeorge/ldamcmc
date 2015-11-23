#' #############################################################################
#' This is a script to test the LDA full Gibbs sampler on the Wikipedia dataset
#' whales  
#' 
#' #############################################################################

# Loads necessary libraries and sets global variables  --------------------

library(ldamcmc)

alpha          <- .1
eta            <- .4
SEED           <- 1983 
max.iter       <- 1000 # the maximum number of Gibbs iterations
burn.in        <- 900
spacing        <- 1
store.z        <- 1                          # store z samples ? 
store.beta     <- 1                          # store beta samples ? 
store.theta    <- 1                          # store theta samples ? 
store.lp       <- 1                          # store log posterior for each iteration 
set.seed(SEED)


# Loads data --------------------------------------------------------------
# Here, the dataset chosen is whales. We can any available data for this purpose 
# e.g., wt, felines, cats, canis, etc. 

data(whales) # See help(whales)

K <- length(class.labels) # the number of topics 
V <- length(vocab) # the vocabulary size 


# Gibbs sampling  ---------------------------------------------------------


alpha.v <- array(alpha, dim=c(K, 1));         

# Full Gibbs sampling (FGS)
# See help(lda_fgs)
# model <- lda_fgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, 
#                  spacing, store.z, store.beta, store.theta, store.lp);

# Augmented Collapsed Gibbs sampling (ACGS)
# NOTE: if store_dirichlet is set as 0, this only do Collapsed Gibbs sampling. 
# See help(lda_acgs)
model <- lda_acgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, 
                  spacing, store.z, store.beta, store.theta, store.lp);


# Displays most probable words from each topic
# See help(calc_top_topic_words)
last.beta.sample <- model$beta[,,max.iter-burn.in]; # takes the last sample 
cat("\nMost probable words from each topic (using the last sample):\n")
calc_top_topic_words(last.beta.sample, vocab, num.words=20, num.digits=2)


# Saves every object into a file  -----------------------------------------


rdata.file <- paste(ds.name, "-gibbs-h(", eta, ",", alpha, ").RData", sep="")[1]
save.image(rdata.file)
cat("\nRData is saved to:", rdata.file)


