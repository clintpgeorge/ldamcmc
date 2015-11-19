#' #############################################################################
#' This is a script to test the LDA full Gibbs sampler on the Wikipedia datasets 
#'
#' Example run:
#'  Rscript lda_gibbs_wikipedia_data.R "../data-raw/" "wt16" "category" .1 .64
#' #############################################################################

## Removes all objects in the current R state 

rm(list=ls());


# Handles commandline arguments  ------------------------------------------

args           <- commandArgs(TRUE)
data.dir       <- args[1] 
prefix         <- args[2] 
class.type     <- args[3] 
alpha          <- as.numeric(args[4]) 
eta            <- as.numeric(args[5])

cat("\nCommandline Arguments:\n")
cat(paste("\ndata.dir:", data.dir))
cat(paste("\nprefix:", prefix))
cat(paste("\nclass.type:", class.type))
cat(paste("\nalpha:", alpha))
cat(paste("\neta:", eta))
cat("\n")



# Loads necessary libraries and sets global variables  --------------------

library(ldamcmc)

# Change the following if you are not passing them thru commandline 

# data.dir       <- "../source/"
# prefix         <- "wt16"
# class.type     <- "category"
# alpha          <- .1
# eta            <- .64

SEED           <- 1983 
max.iter       <- 11000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 1
store.z        <- 1                          # store z samples ? 
store.beta     <- 1                          # store beta samples ? 
store.theta    <- 1                          # store theta samples ? 
store.lp       <- 1                          # store log posterior for each iteration 
set.seed(SEED)


# Loads data --------------------------------------------------------------

cat("\nLoading data files...\n")

setwd(data.dir) # Sets the working dir
vocab.file <- paste(prefix, ".ldac.vocab", sep="")
doc.file <- paste(prefix, ".ldac", sep="")
metadata.file <- paste(prefix, ".csv", sep="")

vocab <- readLines(vocab.file);
documents <- read_docs(doc.file);
doc.metadata <- read.csv2(metadata.file, header=T, sep=';');

ds <- vectorize_docs(documents)
doc.N <- calc_doc_lengths(documents)
num.docs <- nrow(doc.metadata)

class.labels <- levels(doc.metadata[, class.type])
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


rdata.file <- paste(prefix, "-gibbs-h(", eta, ",", alpha, ").RData", sep="")[1]
save.image(rdata.file)
cat("\nRData is saved to:", rdata.file)


