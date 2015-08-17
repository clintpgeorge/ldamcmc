#ifndef LDAMCMC_H
#define LDAMCMC_H

#include <assert.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ;
using namespace std ;

RNGScope scope;

// Helper functions 

/**
 * Samples an integer from [0, K) uniformly at random
 * 
 * Arguments: 
 * 		K - the upper interval 
 * Returns:
 * 		the sampled integer  
 */
unsigned int sample_uniform_int(unsigned int K){
  return (unsigned int) (runif(1)(0) * (double)K); // To speedup
}

/**
* Samples from a given multimomial probability vector
*
* Arguments:
* 		theta - the Multinomial probability vector
* Returns:
* 		t - the sampled index
*
*/
unsigned int sample_multinomial (vec theta) {
  
  unsigned int t = 0;
  double total_prob = accu(theta);
  double u = runif(1)(0) * total_prob;
  double cumulative_prob = theta(0);
  
  while(u > cumulative_prob){
    t++;
    cumulative_prob += theta(t);
  }
  
  return t;
  
}

/**
* Samples from a Dirichlet distribution given a set of hyperparameters
* 
* Aruguments:
* 		num_elements - the dimentionality of the Dirichlet distribution 
* 		alpha - the hyperparameter vector which is in the column vector format  
* Returns: 
* 		the Dirichlet sample in the column vector format   
*/
vec sample_dirichlet(unsigned int num_elements, vec alpha){
  
  vec dirichlet_sample = zeros<vec>(num_elements);
  
  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = R::rgamma(1, alpha(i)); // rgamma(1, alpha(i), 1.0)(0);
  
  dirichlet_sample /= accu(dirichlet_sample);
  
  return dirichlet_sample;
  
}

/**
* Samples from a Dirichlet distribution given a set of hyperparameters
* 
* Aruguments:
* 		num_elements - the dimentionality of the Dirichlet distribution 
* 		alpha - the hyperparameter vector which is in the row vector format  
* Returns: 
* 		the Dirichlet sample in the row vector format   
*/
rowvec sample_dirichlet_row_vec (unsigned int num_elements, rowvec alpha){
  
  rowvec dirichlet_sample = zeros<rowvec>(num_elements);
  
  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = R::rgamma(1, alpha(i)); // rgamma(1, alpha(i), 1.0)(0);
  
  dirichlet_sample /= accu(dirichlet_sample);
  
  return dirichlet_sample;
  
}

/**
* Samples random permutations for a given count
*
* Arguments:
* 		n - the number of samples
* Return:
* 		order - a vector of indices that represents
* 				the permutations of numbers in [1, n]
**/
uvec randperm(unsigned int n) {
  uvec order = zeros<uvec>(n);
  unsigned int k, nn, takeanumber, temp;
  for (k=0; k<n; k++) order(k) = k;
  nn = n;
  for (k=0; k<n; k++) {
    takeanumber = sample_uniform_int(nn); // take a number between 0 and nn-1
    temp = order(nn-1);
    order(nn-1) = order(takeanumber);
    order(takeanumber) = temp;
    nn--;
  }
  return order;
}


vec log_gamma_vec(vec x_vec){
  
  vec lgamma_vec = zeros<vec>(x_vec.n_elem);
  
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    lgamma_vec(i) = lgamma(x_vec(i));
  
  return lgamma_vec;
  
}

vec gamma_col_vec(vec x_vec){
  // It took 2hrs of my time in the April 19, 2014 morning to make this function
  // work. The main issue was with accessing the R gamma function from the 
  // RcppArmadillo namespace. See 
  // http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html 
  // gamma(as<NumericVector>(wrap(x_vec))) is another option, but it seems to be
  // slow. See 
  // http://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector
  
  vec gamma_vec = zeros<vec>(x_vec.n_elem);
  
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    gamma_vec(i) = Rf_gammafn(x_vec(i));
  
  return gamma_vec;
}


/***
* Computes log posterior based on (3.4)
*/
double calc_log_posterior(mat prior_theta, 
                          mat prior_beta,
                          vector < vector < unsigned int > > doc_word_indices, 
                          uvec doc_lengths,
                          uvec word_ids, 
                          uvec z, 
                          vec alpha_v, 
                          double eta){
  
  double lp = 0.0;
  unsigned int d, i, num_docs, num_topics, vocab_size;
  mat log_prior_theta = log(prior_theta);
  mat log_prior_beta = log(prior_beta);
  num_topics = prior_beta.n_rows;
  vocab_size = prior_beta.n_cols;
  num_docs = prior_theta.n_cols;
  
  
  for (d = 0; d < num_docs; d++){ // for each document
    
    vector < unsigned int > word_idx = doc_word_indices[d];
    
    vec n_dj = zeros<vec>(num_topics);
    mat m_djt = zeros<mat>(num_topics, vocab_size);
    for (i = 0; i < doc_lengths(d); i++){
      n_dj(z(word_idx[i])) += 1;
      m_djt(z(word_idx[i]), word_ids(word_idx[i])) += 1;
    }
    
    lp += accu(m_djt % log_prior_beta);
    lp += accu((n_dj + alpha_v - 1.0) % log_prior_theta.col(d));
    
  }
  
  lp += accu((eta - 1.0) * log_prior_beta);
  
  return lp;
}

RcppExport SEXP lda_fgs_lppv(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
                             SEXP doc_lengths_, SEXP topic_assignments_, 
                             SEXP alpha_v_, SEXP eta_, SEXP max_iter_, 
                             SEXP burn_in_, SEXP spacing_);
RcppExport SEXP lda_fgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
                        SEXP doc_lengths_, SEXP topic_assignments_, 
                        SEXP alpha_v_, SEXP eta_, SEXP max_iter_, SEXP burn_in_, 
                        SEXP spacing_, SEXP save_z_, SEXP save_beta_, 
                        SEXP save_theta_, SEXP save_lp_);
RcppExport SEXP lda_fgs_blei_corpus(SEXP num_topics_, SEXP vocab_size_, 
                                    SEXP doc_lengths_, SEXP docs_, 
                                    SEXP topic_assignments_, SEXP alpha_v_, 
                                    SEXP eta_, SEXP max_iter_, SEXP burn_in_, 
                                    SEXP spacing_, SEXP save_z_, 
                                    SEXP save_beta_, SEXP save_theta_, 
                                    SEXP save_lp_);
RcppExport SEXP lda_acgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
                         SEXP doc_lengths_, SEXP topic_assignments_, 
                         SEXP alpha_v_, SEXP eta_, SEXP max_iter_, 
                         SEXP burn_in_, SEXP spacing_, SEXP save_z_, 
                         SEXP save_beta_, SEXP save_theta_, SEXP save_lp_);
  
#endif
