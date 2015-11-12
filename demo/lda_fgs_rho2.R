#' #############################################################################
#' 
#' This script computes estimates of the ratios of the discrepancy rho_2 (Table 
#' 4) in the paper "Principled Selection of Hyperparameters in the Latent 
#' Dirichlet Allocation Model".
#' 
#' Last updated on: November 10, 2015 
#' 
#' Example: 
#' Rscript lda_fgs_rho2.R "D:/workspace/ldamcmc/data-raw" "felines" "category" .5 .155 1983 61000 1000 50
#' Rscript lda_fgs_rho2.R "D:/workspace/ldamcmc/data-raw" "rec" "category" .46 .09 2015 61000 1000 50
#' #############################################################################

## Loads necessary libraries and sets global variables 

rm(list=ls()); # Removes all objects in the current R state 

library(ldamcmc)

# Handles commandline arguments 

args <- commandArgs(TRUE)

data.dir       <- args[1] 
prefix         <- args[2] 
class.type     <- args[3] 
eb.eta         <- as.numeric(args[4])
eb.alpha       <- as.numeric(args[5]) 
SEED           <- as.numeric(args[6]) 
max.iter       <- as.numeric(args[7]) # the maximum number of Gibbs iterations
burn.in        <- as.numeric(args[8]) 
spacing        <- as.numeric(args[9]) 
save.beta      <- 1
save.theta     <- 1

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


## Loads data 

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
K <- length(class.labels)
delta.theta <- array(0.0, c(K, num.docs))
rownames(delta.theta) <- class.labels
for (i in 1:num.docs){
  delta.theta[doc.metadata[i, class.type], i] <- 1.
}

V              <- length(vocab)

rdata.file <- paste(prefix, "-rho2-default-eb-h(", eb.eta, ",", eb.alpha, ")-K", 
                    K, format(Sys.time(), "%Y%b%d%H%M%S"), ".RData", sep="")[1]


## default hyperparameters 

# gensim package (2010)
d1.alpha.v <- array(1/K, dim=c(K, 1));                          
d1.eta <- 1/K;
set.seed(SEED)
d1.model <- lda_fgs(K, V, ds$wid+1, doc.N, d1.alpha.v, d1.eta, max.iter, 
                    burn.in, spacing, save.theta=save.theta, save.beta=save.beta);

# Asuncion et al. (2009) 
d2.alpha.v <- array(0.1, dim=c(K, 1));                          
d2.eta <- 0.1;
set.seed(SEED)
d2.model <- lda_fgs(K, V, ds$wid+1, doc.N, d2.alpha.v, d2.eta, max.iter, 
                    burn.in, spacing, save.theta=save.theta, save.beta=save.beta);

# Griffiths and Steyvers (2004)
d3.alpha.v <- array(50/K, dim=c(K, 1));                          
d3.eta <- 0.1;
set.seed(SEED)
d3.model <- lda_fgs(K, V, ds$wid+1, doc.N, d3.alpha.v, d3.eta, max.iter, 
                    burn.in, spacing, save.theta=save.theta, save.beta=save.beta);

## Based on empirical Bayes approach 
eb.alpha.v <- array(eb.alpha, dim=c(K, 1)); 
set.seed(SEED)
eb.model <- lda_fgs(K, V, ds$wid+1, doc.N, eb.alpha.v, eb.eta, max.iter, 
                    burn.in, spacing, save.theta=save.theta, save.beta=save.beta);


## Saves every object into a file 
save.image(rdata.file)

# Measure distances between document's true topic and estimated topic 
# distributions 

ctf <- calc_class_term_frequency(class.labels, vocab, 
                                 doc.metadata[, class.type], documents);
nctf <- normalize(ctf, dim=2);
num.samples <- dim(eb.model$theta)[3];

eb.topic.labels <- calc_beta_topic_labels(nctf, eb.model$beta)
d1.topic.labels <- calc_beta_topic_labels(nctf, d1.model$beta)
d2.topic.labels <- calc_beta_topic_labels(nctf, d2.model$beta)
d3.topic.labels <- calc_beta_topic_labels(nctf, d3.model$beta)


# with grouping 

btheta.L1 <- array(0., c(num.samples, 4));
dim.names <- list(rownames(delta.theta), 1:num.docs);
for (i in 1:num.samples){
  eb.theta <- array(0., c(K, num.docs), dimnames=dim.names);
  d1.theta <- array(0., c(K, num.docs), dimnames=dim.names);
  d2.theta <- array(0., c(K, num.docs), dimnames=dim.names);
  d3.theta <- array(0., c(K, num.docs), dimnames=dim.names);
  for (j in 1:K){
    eb.theta[eb.topic.labels[j,i],] <- eb.theta[eb.topic.labels[j,i],] + eb.model$theta[j,,i];   
    d1.theta[d1.topic.labels[j,i],] <- d1.theta[d1.topic.labels[j,i],] + d1.model$theta[j,,i];
    d2.theta[d2.topic.labels[j,i],] <- d2.theta[d2.topic.labels[j,i],] + d2.model$theta[j,,i];
    d3.theta[d3.topic.labels[j,i],] <- d3.theta[d3.topic.labels[j,i],] + d3.model$theta[j,,i];
  }
  btheta.L1[i,1] <- sum(abs(eb.theta - delta.theta));
  btheta.L1[i,2] <- sum(abs(d1.theta - delta.theta));
  btheta.L1[i,3] <- sum(abs(d2.theta - delta.theta));
  btheta.L1[i,4] <- sum(abs(d3.theta - delta.theta));
}

eb.h.title <- paste("eb-h(", eb.eta, ",", eb.alpha, ")", sep="")[1];
colnames(btheta.L1) <- c(eb.h.title, 'h=(1/K,1/K)', 'h=(.1,.1)', 'h=(.1,50/K)');
rho.2 <- colMeans(btheta.L1) # rho.2 (4.2)
rho.2
rho.2/rho.2[1]



## Saves every object into a file 

save.image(rdata.file)
cat("\nRData is saved to:", rdata.file, "\n")




