#ifndef LDAMCMC_H
#define LDAMCMC_H

#include <assert.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ;
using namespace std ;

RNGScope scope;

// Helper functions 

unsigned int sample_multinomial (vec theta);
vec sample_dirichlet (unsigned int num_elements, vec alpha);
rowvec sample_dirichlet_row_vec (unsigned int num_elements, rowvec alpha);
double calc_log_posterior(mat prior_theta, mat prior_beta,
  vector < vector < unsigned int > > doc_word_indices, uvec doc_lengths,
	uvec word_ids, uvec z, vec alpha_v, double eta);


// Main functions 

RcppExport SEXP lda_fgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
  SEXP doc_lengths_, SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_, 
  SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_);
  
RcppExport SEXP lda_fgs_blei_corpus(SEXP num_topics_, SEXP vocab_size_, 
  SEXP doc_lengths_, SEXP docs_, SEXP topic_assignments_, SEXP alpha_v_, 
  SEXP eta_, SEXP max_iter_, SEXP burn_in_, SEXP spacing_, 
  SEXP store_dirichlet_);

RcppExport SEXP lda_acgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
  SEXP doc_lengths_, SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_, 
  SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_);


#endif
