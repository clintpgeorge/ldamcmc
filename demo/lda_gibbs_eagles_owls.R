#' #############################################################################
#' This is a script to test the LDA full Gibbs sampler on the Wikipedia dataset
#' eagles-owls  
#' 
#' #############################################################################

# Loads necessary libraries and sets global variables  --------------------

library(ldamcmc)

alpha          <- .1
eta            <- .4
SEED           <- 1983 
max.iter       <- 5100 # the maximum number of Gibbs iterations
burn.in        <- 5000
spacing        <- 1
store.z        <- 1                          # store z samples ? 
store.beta     <- 1                          # store beta samples ? 
store.theta    <- 1                          # store theta samples ? 
store.lp       <- 1                          # store log posterior for each iteration 
set.seed(SEED)


# Loads data --------------------------------------------------------------
# Here, the dataset chosen is eagles-owls. We can any available data for this 
# purpose 
# e.g., wt, felines, cats, canis, etc. 

data("eagles-owls") # See help("eagles-owls")

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


last.theta.sample <- model$theta[,,max.iter-burn.in]; # takes the last sample 
x <- last.theta.sample[1,]
y <- last.theta.sample[2,]
color.map <- rep("red", num.docs)
color.map[docs.metadata[,"category"] == "Category:Eagles"] <- "blue"

# Sets cex properties for the line plots 
cex.axis <- 1.3;
cex.lab <- 1.5;
cex <- 1; 
op <- par(bg = "white")
trellis.device(pdf, file="eagles-owls-lda-space.pdf", height=5.5, width=7, onefile=T)
plot(x,y, main="LDA Topic Space", xlab="Topic-1", ylab="Topic-2", 
     col=color.map, pch=16, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis);
par(op)
dev.off() 


# Saves every object into a file  -----------------------------------------


rdata.file <- paste(ds.name, "-gibbs-h(", eta, ",", alpha, ").RData", sep="")[1]
save.image(rdata.file)
cat("\nRData is saved to:", rdata.file)

