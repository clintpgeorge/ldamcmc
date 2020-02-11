// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lda_acgs_st
List lda_acgs_st(unsigned int num_topics, unsigned int vocab_size, List docs_tf, arma::mat h_grid, arma::mat st_grid, List st_grid_nbrs, unsigned int init_st_grid_index, arma::vec zetas, unsigned int tuning_iter, unsigned int max_iter_tuning, unsigned int max_iter_final, unsigned int burn_in, unsigned int spacing, double test_set_share, bool save_beta, bool save_theta, bool save_lp, bool save_hat_ratios, bool save_tilde_ratios, int verbose);
RcppExport SEXP _ldamcmc_lda_acgs_st(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP h_gridSEXP, SEXP st_gridSEXP, SEXP st_grid_nbrsSEXP, SEXP init_st_grid_indexSEXP, SEXP zetasSEXP, SEXP tuning_iterSEXP, SEXP max_iter_tuningSEXP, SEXP max_iter_finalSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP test_set_shareSEXP, SEXP save_betaSEXP, SEXP save_thetaSEXP, SEXP save_lpSEXP, SEXP save_hat_ratiosSEXP, SEXP save_tilde_ratiosSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_grid(h_gridSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type st_grid(st_gridSEXP);
    Rcpp::traits::input_parameter< List >::type st_grid_nbrs(st_grid_nbrsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type init_st_grid_index(init_st_grid_indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type zetas(zetasSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type tuning_iter(tuning_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter_tuning(max_iter_tuningSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter_final(max_iter_finalSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< bool >::type save_hat_ratios(save_hat_ratiosSEXP);
    Rcpp::traits::input_parameter< bool >::type save_tilde_ratios(save_tilde_ratiosSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_acgs_st(num_topics, vocab_size, docs_tf, h_grid, st_grid, st_grid_nbrs, init_st_grid_index, zetas, tuning_iter, max_iter_tuning, max_iter_final, burn_in, spacing, test_set_share, save_beta, save_theta, save_lp, save_hat_ratios, save_tilde_ratios, verbose));
    return rcpp_result_gen;
END_RCPP
}
// lda_cgs_perplexity
List lda_cgs_perplexity(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_theta, bool save_beta, bool save_lp, int verbose, double test_set_share);
RcppExport SEXP _ldamcmc_lda_cgs_perplexity(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP test_set_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_cgs_perplexity(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in, spacing, save_theta, save_beta, save_lp, verbose, test_set_share));
    return rcpp_result_gen;
END_RCPP
}
// lda_cgs_em_perplexity
List lda_cgs_em_perplexity(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned int em_max_iter, unsigned int gibbs_max_iter, unsigned int burn_in, unsigned int spacing, bool save_theta, bool save_beta, bool save_lp, int verbose, double test_set_share);
RcppExport SEXP _ldamcmc_lda_cgs_em_perplexity(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP em_max_iterSEXP, SEXP gibbs_max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP test_set_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type em_max_iter(em_max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type gibbs_max_iter(gibbs_max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_cgs_em_perplexity(num_topics, vocab_size, docs_tf, alpha_h, eta_h, em_max_iter, gibbs_max_iter, burn_in, spacing, save_theta, save_beta, save_lp, verbose, test_set_share));
    return rcpp_result_gen;
END_RCPP
}
// lda_cgs_em
List lda_cgs_em(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned long int em_max_iter, unsigned long int gibbs_max_iter, unsigned long int burn_in, unsigned long int spacing, int verbose);
RcppExport SEXP _ldamcmc_lda_cgs_em(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP em_max_iterSEXP, SEXP gibbs_max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned long int >::type em_max_iter(em_max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned long int >::type gibbs_max_iter(gibbs_max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned long int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned long int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_cgs_em(num_topics, vocab_size, docs_tf, alpha_h, eta_h, em_max_iter, gibbs_max_iter, burn_in, spacing, verbose));
    return rcpp_result_gen;
END_RCPP
}
// lda_fgs_perplexity
List lda_fgs_perplexity(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_theta, bool save_beta, bool save_lp, int verbose, double test_set_share);
RcppExport SEXP _ldamcmc_lda_fgs_perplexity(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP verboseSEXP, SEXP test_set_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_fgs_perplexity(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in, spacing, save_theta, save_beta, save_lp, verbose, test_set_share));
    return rcpp_result_gen;
END_RCPP
}
// lda_fgs_BF_perplexity
List lda_fgs_BF_perplexity(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, arma::mat h_grid, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, bool save_theta, bool save_beta, bool save_lp, bool save_BF, int verbose, double test_set_share);
RcppExport SEXP _ldamcmc_lda_fgs_BF_perplexity(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP h_gridSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP save_thetaSEXP, SEXP save_betaSEXP, SEXP save_lpSEXP, SEXP save_BFSEXP, SEXP verboseSEXP, SEXP test_set_shareSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_grid(h_gridSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< bool >::type save_BF(save_BFSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_fgs_BF_perplexity(num_topics, vocab_size, docs_tf, alpha_h, eta_h, h_grid, max_iter, burn_in, spacing, save_theta, save_beta, save_lp, save_BF, verbose, test_set_share));
    return rcpp_result_gen;
END_RCPP
}
// lda_fgs_ppc
List lda_fgs_ppc(unsigned int num_topics, unsigned int vocab_size, List docs_tf, double alpha_h, double eta_h, unsigned int max_iter, unsigned int burn_in, unsigned int spacing, int verbose);
RcppExport SEXP _ldamcmc_lda_fgs_ppc(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP alpha_hSEXP, SEXP eta_hSEXP, SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_h(alpha_hSEXP);
    Rcpp::traits::input_parameter< double >::type eta_h(eta_hSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_fgs_ppc(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in, spacing, verbose));
    return rcpp_result_gen;
END_RCPP
}
// lda_fgs_st_perplexity
List lda_fgs_st_perplexity(unsigned int num_topics, unsigned int vocab_size, List docs_tf, arma::mat h_grid, arma::mat st_grid, List st_grid_nbrs, unsigned int init_st_grid_index, arma::vec zetas, unsigned int tuning_iter, unsigned int max_iter_tuning, unsigned int max_iter_final, unsigned int burn_in, unsigned int spacing, double test_set_share, bool save_beta, bool save_theta, bool save_lp, bool save_hat_ratios, bool save_tilde_ratios, int verbose);
RcppExport SEXP _ldamcmc_lda_fgs_st_perplexity(SEXP num_topicsSEXP, SEXP vocab_sizeSEXP, SEXP docs_tfSEXP, SEXP h_gridSEXP, SEXP st_gridSEXP, SEXP st_grid_nbrsSEXP, SEXP init_st_grid_indexSEXP, SEXP zetasSEXP, SEXP tuning_iterSEXP, SEXP max_iter_tuningSEXP, SEXP max_iter_finalSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP test_set_shareSEXP, SEXP save_betaSEXP, SEXP save_thetaSEXP, SEXP save_lpSEXP, SEXP save_hat_ratiosSEXP, SEXP save_tilde_ratiosSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_topics(num_topicsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type vocab_size(vocab_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type docs_tf(docs_tfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type h_grid(h_gridSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type st_grid(st_gridSEXP);
    Rcpp::traits::input_parameter< List >::type st_grid_nbrs(st_grid_nbrsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type init_st_grid_index(init_st_grid_indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type zetas(zetasSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type tuning_iter(tuning_iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter_tuning(max_iter_tuningSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter_final(max_iter_finalSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< double >::type test_set_share(test_set_shareSEXP);
    Rcpp::traits::input_parameter< bool >::type save_beta(save_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_theta(save_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type save_lp(save_lpSEXP);
    Rcpp::traits::input_parameter< bool >::type save_hat_ratios(save_hat_ratiosSEXP);
    Rcpp::traits::input_parameter< bool >::type save_tilde_ratios(save_tilde_ratiosSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_fgs_st_perplexity(num_topics, vocab_size, docs_tf, h_grid, st_grid, st_grid_nbrs, init_st_grid_index, zetas, tuning_iter, max_iter_tuning, max_iter_final, burn_in, spacing, test_set_share, save_beta, save_theta, save_lp, save_hat_ratios, save_tilde_ratios, verbose));
    return rcpp_result_gen;
END_RCPP
}
// sample_antoniak
double sample_antoniak(unsigned int N, double alpha);
RcppExport SEXP _ldamcmc_sample_antoniak(SEXP NSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_antoniak(N, alpha));
    return rcpp_result_gen;
END_RCPP
}
// sample_multinomial
unsigned int sample_multinomial(arma::vec theta);
RcppExport SEXP _ldamcmc_sample_multinomial(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_multinomial(theta));
    return rcpp_result_gen;
END_RCPP
}
// sample_dirichlet
arma::vec sample_dirichlet(unsigned int num_elements, arma::vec alpha);
RcppExport SEXP _ldamcmc_sample_dirichlet(SEXP num_elementsSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type num_elements(num_elementsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_dirichlet(num_elements, alpha));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP lda_acgs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP lda_acgs_hs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP lda_fgs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP lda_fgs_hs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP lda_fgs_lppv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP lda_fgs_st(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ldamcmc_lda_acgs_st", (DL_FUNC) &_ldamcmc_lda_acgs_st, 20},
    {"_ldamcmc_lda_cgs_perplexity", (DL_FUNC) &_ldamcmc_lda_cgs_perplexity, 13},
    {"_ldamcmc_lda_cgs_em_perplexity", (DL_FUNC) &_ldamcmc_lda_cgs_em_perplexity, 14},
    {"_ldamcmc_lda_cgs_em", (DL_FUNC) &_ldamcmc_lda_cgs_em, 10},
    {"_ldamcmc_lda_fgs_perplexity", (DL_FUNC) &_ldamcmc_lda_fgs_perplexity, 13},
    {"_ldamcmc_lda_fgs_BF_perplexity", (DL_FUNC) &_ldamcmc_lda_fgs_BF_perplexity, 15},
    {"_ldamcmc_lda_fgs_ppc", (DL_FUNC) &_ldamcmc_lda_fgs_ppc, 9},
    {"_ldamcmc_lda_fgs_st_perplexity", (DL_FUNC) &_ldamcmc_lda_fgs_st_perplexity, 20},
    {"_ldamcmc_sample_antoniak", (DL_FUNC) &_ldamcmc_sample_antoniak, 2},
    {"_ldamcmc_sample_multinomial", (DL_FUNC) &_ldamcmc_sample_multinomial, 1},
    {"_ldamcmc_sample_dirichlet", (DL_FUNC) &_ldamcmc_sample_dirichlet, 2},
    {"lda_acgs",     (DL_FUNC) &lda_acgs,     14},
    {"lda_acgs_hs",  (DL_FUNC) &lda_acgs_hs,  16},
    {"lda_fgs",      (DL_FUNC) &lda_fgs,      14},
    {"lda_fgs_hs",   (DL_FUNC) &lda_fgs_hs,   16},
    {"lda_fgs_lppv", (DL_FUNC) &lda_fgs_lppv, 10},
    {"lda_fgs_st",   (DL_FUNC) &lda_fgs_st,   23},
    {NULL, NULL, 0}
};

RcppExport void R_init_ldamcmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
