# include "utils.h"

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

extern
double lda_log_posterior(
    arma::uvec doc_word_counts,
    arma::mat theta_samples,
    arma::mat beta_samples,
    vector < vector < unsigned int > > doc_word_ids,
    vector < vector < unsigned int > > doc_word_zids,
    double alpha_h,
    double eta_h
);

extern 
double calc_log_posterior(
    mat prior_theta,
    mat prior_beta,
    vector < vector < unsigned int > > doc_word_indices,
    uvec doc_lengths,
    uvec word_ids,
    uvec z,
    vec alpha_v,
    double eta
); 

extern 
double calc_prior_ratio(
    mat beta, 
    mat theta, 
    double alpha, 
    double eta, 
    double base_alpha, 
    double base_eta
); 

