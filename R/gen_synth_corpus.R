calc_topic_counts <- function(Z, K)
{
  Nt <- array(0, c(1, K)); 
  for (k in 1:K) Nt[k] <- sum(Z == k);
  
  return(Nt); 	
}

#' Generates a synthetic corpus 
#' 
#' Generates document words using the LDA generative process given a beta. It's
#' used to test the correctness of the Gibbs sampling algorithms.
#' 
#' @param D the number of documents in the corpus
#' @param lambda.hat the mean of document counts
#' @param alpha.v the vector of Dirichlet hyperparameters (K X 1) for document
#'   topic mixtures
#' @param beta the beta matrix (counts) for topic word probabilities (K x V
#'   format)
#'   
#' @return a list of generated documents' details
#'   
#' @export
#' 
#' @examples
#' K              <- 2
#' V              <- 20 
#' D              <- 100
#' gen.alpha.v    <- array(7, c(K, 1)); 
#' gen.eta.v      <- array(3, c(1, V)); 
#' lambda.hat     <- 80
#' 
#' ## Generates the synthetic beta.m
#' beta.m         <- matrix(1e-2, nrow=K, ncol=V)
#' beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
#' beta.m[2, ]    <- rdirichlet(1, gen.eta.v);
#' 
#' ## Generates documents with a given beta.m
#' ds             <- gen_synth_corpus(D, lambda.hat, gen.alpha.v, beta.m);
#' 
gen_synth_corpus <- function(D, lambda.hat, alpha.v, beta)
{
  K <- nrow(beta);                                 # the number of topics 
  V <- ncol(beta)                                  # the vocabulary size 
  theta.counts <- matrix(0, nrow=K, ncol=D);       # stores document topic word counts  
  beta.counts <- matrix(0, nrow=K, ncol=V);        # stores topic word counts 
  theta.samples <- matrix(0, nrow=K, ncol=D);      
  
  did <- c();
  wid <- c();
  zid <- c();
  doc.N <- array(lambda.hat, dim=c(D, 1)); 
  
  num <- 1;
  doc.idx <- vector("list", D);
  
  for (d in 1:D)
  {
    ptm <- proc.time();
    
    theta.samples[,d] <- rdirichlet(1, alpha.v);
    
    did <- cbind(did, array(1, c(1, doc.N[d])) * d); # document instances 
    
    z_d <- c(); 
    indices <- c();
    for (i in 1:doc.N[d]){
      
      z_dn <- which(rmultinom(1, size=1, prob=theta.samples[,d]) == 1); # samples topic 
      w_dn <- which(rmultinom(1, size=1, beta[z_dn,]) == 1); # samples word
      
      wid <- cbind(wid, w_dn); 
      z_d <- cbind(z_d, z_dn);  
      indices <- cbind(indices, num);
      
      num <- num + 1;            
      
    }
    doc.idx[[d]] <- as.integer(indices); # stores the document word indices     
    
    
    theta.counts[, d] <- calc_topic_counts(z_d, K); # calculates the document topic counts
    
    zid <- cbind(zid, z_d); 
    
    ptm <- proc.time() - ptm;
    cat("document = ", d, " time = ", ptm[3], " # words = ", doc.N[d], "\n");
  }
  
  
  total.N <- sum(doc.N);
  
  for (i in 1:total.N){ 
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1;
  }
  
  list(did=as.vector(did), wid=as.vector(wid), zid=as.vector(zid), 
    theta.counts=theta.counts, beta.counts=beta.counts, 
    theta.samples=theta.samples, total.N=total.N, doc.N=doc.N, 
    doc.idx=doc.idx);
  
}

#' Generates a synthetic corpus based on symmetric Dirichlets 
#' 
#' Generates documents using the LDA generative process based on a set of 
#' predefined values. 
#' 
#' @param K number of topics
#' @param V vocabulary size 
#' @param D number of documents
#' @param doc.size number of words in each document
#' @param alpha hyperparameter for document Dirichlet sampling 
#' @param eta hyperparameter for topic Diriclet sampling 
#'   
#' @return a list of generated documents and their statistics 
#'   
#' @export
#' 
#' @family corpus 
#' 
#' @details Last modified on: May 24, 2015 
#' 
#' @examples
#
#' ## Generates documents with given parameters 
#' 
#' K              <- 2
#' V              <- 20 
#' D              <- 100
#' alpha          <- 7
#' eta            <- 3
#' doc.size       <- 80
#' 
#' ds             <- gen_corpus(K, V, D, doc.size, alpha, eta);
#' 
gen_corpus <- function(K, V, D, doc.size, alpha, eta)
{
  theta.counts <- matrix(0, nrow=K, ncol=D);       # stores document topic word counts  
  beta.counts <- matrix(0, nrow=K, ncol=V);        # stores topic word counts 
  theta.samples <- matrix(0, nrow=K, ncol=D);      
  beta.samples <- matrix(1e-2, nrow=K, ncol=V); 
  
  alpha.v <- array(alpha, c(K, 1)); 
  eta.v <- array(eta, c(1, V)); 
  did <- c();
  wid <- c();
  zid <- c();
  doc.N <- array(doc.size, dim=c(D, 1)); 
  doc.idx <- vector("list", D);
  word.idx <- 1; # initialize the corpus word index
  
  ptm <- proc.time();
  
  for (j in 1:K){ # Topic Dirichlet sampling 
    beta.samples[j, ]  <- rdirichlet(1, eta.v);
  }
  
  for (d in 1:D){ # Document sampling 
    
    theta.samples[,d] <- rdirichlet(1, alpha.v);
    
    did <- cbind(did, array(1, c(1, doc.N[d])) * d); # document instances 
    z_d <- c(); 
    indices <- c();
    
    for (i in 1:doc.N[d]){ # Word sampling 
      
      z_dn <- which(rmultinom(1, size=1, prob=theta.samples[,d]) == 1); # samples topic 
      w_dn <- which(rmultinom(1, size=1, beta.samples[z_dn,]) == 1); # samples word
      
      wid <- cbind(wid, w_dn); 
      z_d <- cbind(z_d, z_dn);  
      indices <- cbind(indices, word.idx);
      word.idx <- word.idx + 1;            
      
    }
    
    doc.idx[[d]] <- as.integer(indices); # stores the document word indices     
    theta.counts[, d] <- calc_topic_counts(z_d, K); # calculates the document topic counts
    zid <- cbind(zid, z_d); 
    
  }
  
  total.N <- sum(doc.N);
  for (i in 1:total.N){ 
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1;
  }
  
  ptm <- proc.time() - ptm;
  cat("Corpus generation time: ", ptm[3], ", number of total words: ", 
      total.N, "\n", sep="");  
  
  # returns a list 
  list(did=as.vector(did), 
       wid=as.vector(wid), 
       zid=as.vector(zid), 
       theta.counts=theta.counts, 
       beta.counts=beta.counts, 
       theta.samples=theta.samples, 
       beta.samples=beta.samples,
       total.N=total.N, 
       doc.N=doc.N, 
       doc.idx=doc.idx);
}

