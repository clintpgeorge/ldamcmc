#' #############################################################################
#' This is a script to test the LDA full Gibbs sampler on the Wikipedia dataset
#' whales  
#'
#' #############################################################################

# Loads necessary libraries and sets global variables  --------------------

rm(list=ls()); # Removes all objects in the current R state 
library(ldamcmc)


prefix         <- "whales"
class.type     <- "category"
alpha          <- .1
eta            <- .4
SEED           <- 1983 
max.iter       <- 1000 # the maximum number of Gibbs iterations
burn.in        <- 900
spacing        <- 1
store.Dir      <- 1
set.seed(SEED)


# Loads data --------------------------------------------------------------


# See help(whales)
data(whales.vocab)
data(whales.docs)
data(whales.docs.metadata)

ds <- vectorize_docs(whales.docs)
doc.N <- calc_doc_lengths(whales.docs)
num.docs <- nrow(whales.docs.metadata)

class.labels <- levels(whales.docs.metadata[, class.type])
K <- length(class.labels) # the number of topics 
V <- length(whales.vocab) # the vocabulary size 



# Gibbs sampling  ---------------------------------------------------------


alpha.v <- array(alpha, dim=c(K, 1));         

# Full Gibbs sampling (FGS)
# See help(lda_fgs)
model <- lda_fgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, 
                 spacing, store.Dir);

# Augmented Collapsed Gibbs sampling (ACGS)
# NOTE: if store_dirichlet is set as 0, this only do Collapsed Gibbs sampling. 
# See help(lda_acgs)
# model <- lda_acgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, 
#                   spacing, store.Dir);


# Displays most probable words from each topic
# See help(calc_top_topic_words)
last.beta.sample <- model$beta[,,max.iter-burn.in]; # takes the last sample 
cat("\nMost probable words from each topic (using the last sample):\n")
calc_top_topic_words(last.beta.sample, whales.vocab, num.words=20, num.digits=2)


# Saves every object into a file  -----------------------------------------


rdata.file <- paste(prefix, "-gibbs-h(", eta, ",", alpha, ").RData", sep="")[1]
save.image(rdata.file)
cat("\nRData is saved to:", rdata.file)


