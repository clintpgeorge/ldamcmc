// Deprecated on: May 16, 2017 
// 
// Note: All functions in this file use a common, less efficient, format for the 
// input corpus. All these functions are kept to support backward compatibility
//  
#include "utils.h"
#include "lda.h"

#ifndef LDAMCMC_H
#define LDAMCMC_H

RcppExport SEXP lda_fgs(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP word_ids_, 
    SEXP doc_lengths_, 
    SEXP topic_assignments_, 
    SEXP alpha_v_, 
    SEXP eta_, 
    SEXP max_iter_, 
    SEXP burn_in_, 
    SEXP spacing_, 
    SEXP save_z_, 
    SEXP save_beta_, 
    SEXP save_theta_, 
    SEXP save_lp_
);

RcppExport SEXP lda_fgs_lppv(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP word_ids_, 
    SEXP doc_lengths_, 
    SEXP topic_assignments_,
    SEXP alpha_v_,
    SEXP eta_,
    SEXP max_iter_, 
    SEXP burn_in_,
    SEXP spacing_
);

 
 
RcppExport SEXP lda_acgs(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP word_ids_, 
    SEXP doc_lengths_, 
    SEXP topic_assignments_, 
    SEXP alpha_v_, 
    SEXP eta_, 
    SEXP max_iter_, 
    SEXP burn_in_, 
    SEXP spacing_, 
    SEXP save_z_, 
    SEXP save_beta_,
    SEXP save_theta_, 
    SEXP save_lp_
);


RcppExport SEXP lda_fgs_st(
  SEXP num_topics_, 
  SEXP vocab_size_, 
  SEXP word_ids_, 
  SEXP doc_lengths_, 
  SEXP topic_assignments_, 
  SEXP h_grid_, 
  SEXP st_grid_, 
  SEXP st_grid_nbrs_, 
  SEXP init_st_grid_index_, 
  SEXP init_st_grid_zetas_, 
  SEXP max_iter_, 
  SEXP burn_in_, 
  SEXP spacing_, 
  SEXP tuning_iter_, 
  SEXP save_z_, 
  SEXP save_beta_, 
  SEXP save_theta_, 
  SEXP save_st_grid_index_, 
  SEXP save_lp_,
  SEXP save_hat_ratios_,
  SEXP save_tilde_ratios_,
  SEXP verbose_,
  SEXP max_iter_final_
  );

RcppExport SEXP lda_fgs_hs(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP word_ids_, 
    SEXP doc_lengths_, 
    SEXP topic_assignments_, 
    SEXP alpha_, 
    SEXP eta_, 
    SEXP h_grid_, 
    SEXP max_iter_, 
    SEXP burn_in_, 
    SEXP spacing_, 
    SEXP save_z_, 
    SEXP save_beta_, 
    SEXP save_theta_, 
    SEXP save_Bh_, 
    SEXP save_lp_
);

RcppExport SEXP lda_acgs_hs(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP word_ids_, 
    SEXP doc_lengths_, 
    SEXP topic_assignments_, 
    SEXP alpha_, 
    SEXP eta_, 
    SEXP h_grid_, 
    SEXP max_iter_, 
    SEXP burn_in_, 
    SEXP spacing_, 
    SEXP save_z_, 
    SEXP save_beta_, 
    SEXP save_theta_, 
    SEXP save_Bh_, 
    SEXP save_lp_
);
						 
#endif
