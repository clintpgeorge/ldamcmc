# include "lda.h"

///////////////////////////////////////////////////////////////////////////////
// LDA: Helper functions
///////////////////////////////////////////////////////////////////////////////

/***
* Computes the log posterior probability of the LDA model up to a
* multiplicative constant.
*
*/
double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    double alpha_h,
    double eta_h){

  double lp = 0.0;
  unsigned int d, i, num_docs, num_topics, vocab_size;
  arma::mat log_theta = log(theta_samples);
  arma::mat log_beta = log(beta_samples);
  num_topics = beta_samples.n_rows;
  vocab_size = beta_samples.n_cols;
  num_docs = theta_samples.n_cols;
  arma::vec n_dj(num_topics);
  arma::mat m_djt(num_topics, vocab_size);

  for (d = 0; d < num_docs; d++){ // for each document

    vector < unsigned int > word_ids = doc_word_ids[d];
    vector < unsigned int > word_zids = doc_word_zids[d];

    n_dj.fill(0.);
    m_djt.fill(0.);
    for (i = 0; i < doc_word_counts(d); i++){
      n_dj(word_zids[i]) += 1;
      m_djt(word_zids[i], word_ids[i]) += 1;
    }

    lp += arma::accu(m_djt % log_beta);
    lp += arma::accu((n_dj + alpha_h - 1.0) % log_theta.col(d));

  }

  lp += arma::accu((eta_h - 1.0) * log_beta);

  return lp;
}


/***
 * Computes log posterior based on (3.4)
 */
double calc_log_posterior(
    mat prior_theta,
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


/**
* Computes \nu_h(\psi)/\nu_{h_*}(\psi) for a given h, h_*, and 
* \psi = (\beta, \theta, z)
* 
* Reference: 
*  P. George and Doss (2015), Hyperparameter Selection in the Latent Dirichlet 
*  Allocation Model, Equation 2.3 
*/
double calc_prior_ratio(
    mat beta, 
    mat theta, 
    double alpha, 
    double eta, 
    double base_alpha, 
    double base_eta
){
  double num_docs = theta.n_cols; 
  double num_topics = theta.n_rows;
  double vocab_size = beta.n_cols; 
  double ln_dir_ratio = num_docs * (lgamma(num_topics * alpha) 
    - (num_topics * lgamma(alpha)) 
    + (num_topics * lgamma(base_alpha)) 
    - lgamma(num_topics * base_alpha)) + 
    num_topics * (lgamma(vocab_size * eta) 
                    - (vocab_size * lgamma(eta)) 
                    + (vocab_size * lgamma(base_eta)) 
                    - lgamma(vocab_size * base_eta)); 
                    double ln_bf = accu((alpha - base_alpha) * log(theta)) 
                      + accu((eta - base_eta) * log(beta)) 
                      + ln_dir_ratio; 
                      return exp(ln_bf);
}






