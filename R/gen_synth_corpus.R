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
