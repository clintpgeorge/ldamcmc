// =============================================================================
// Serial Tempering 
// =============================================================================

# include "utils.h"
# include "lda.h"

/**
 * Extend division reminder to vectors
 *
 * @param   a       Dividend
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n){
  return a - floor(a/n)*n;
}


//' LDA: Serial Tempering with Perplexity Computation
//'
//' Implements the LDA serial tempering algorithm. Sampling \code{z_{di}}'s 
//' is adapted from the idea of collapsed Gibbs sampling chain (Griffiths and 
//' Steyvers, 2004). To compute perplexity, it first partitions each document in 
//' the corpus into two sets of words: 
//'   (a) a test set (held-out set) and 
//'   (b) a training set, given a user defined \eqn{test_set_share}. 
//' Then, it runs the Markov chain based on the training set and computes 
//' perplexity for the held-out set.
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_tf A list of corpus documents read from the Blei corpus using 
//'   \code{\link{read_docs}} (term indices starts with 0)
//' @param h_grid A 2-dimensional grid of hyperparameters \eqn{h = (\eta, 
//'   \alpha)}. It is a 2 x G matrix, where G is the number of grid points and 
//'   the first row is for \eqn{\alpha} values and the second row is for 
//'   \eqn{\eta} values
//' @param st_grid A 2-dimensional grid of hyperparameters \eqn{h = (\eta, 
//'   \alpha)}. It is a 2 x G matrix, where G is the number of grid points and 
//'   the first row is for \eqn{\alpha} values and the second row is for 
//'   \eqn{\eta} values. This a subgrid on h_grid_ that is used for Serial 
//'   Tempering
//' @param st_grid_nbrs The neighbor indices, from [0, G-1], of each helper grid
//'   point
//' @param init_st_grid_index Index of the helper h grid, from [1, G], of the 
//'   initial hyperparameter \eqn{h = (\eta, \alpha)}
//' @param zetas  Initial guess for normalization constants
//' @param tuning_iter Number of tuning iterations
//' @param max_iter_tuning Maximum number of Gibbs iterations to be performed
//'   for the tuning iterations
//' @param max_iter_final Maximum number of Gibbs iterations to be performed for
//'   the final run
//' @param burn_in Burn-in-period for the Gibbs sampler
//' @param spacing Spacing between the stored samples (to reduce correlation)
//' @param test_set_share Proportion of the test words in each document. Must be
//'   between 0. and 1.
//' @param save_beta If 0 the function does not save \eqn{\beta} samples
//' @param save_theta If 0 the function does not save \eqn{\theta} samples
//' @param save_lp if 0 The function does not save computed log posterior for 
//'   iterations
//' @param save_hat_ratios If 0 the function does not save hat ratios for 
//'   iterations
//' @param save_tilde_ratios If 0 the function does not save tilde ratios for 
//'   iterations
//' @param verbose Values from {0, 1, 2}
//'
//' @return A list of
//'   \item{corpus_topic_counts}{corpus-level topic counts from last iteration
//'   of the Markov chain}
//'   \item{theta_counts}{document-level topic counts from last iteration
//'   of the Markov chain}
//'   \item{beta_counts}{topic word counts from last iteration of the Markov chain}
//'   \item{theta_samples}{\eqn{\theta} samples after the burn in period, if
//'   \code{save_theta} is set}
//'   \item{beta_samples}{\eqn{\beta} samples after the burn in period, if
//'   \code{save_beta} is set}
//'   \item{log_posterior}{the log posterior (upto a constant multiplier) of
//'   the hidden variable \eqn{\psi = (\beta, \theta, z)} in the LDA model,
//'   if \code{save_lp} is set}
//'   \item{perplexity}{perplexity of the held-out words' set}
//'
//' @export
//'
//' @family MCMC
//' 
//' @note 
//'  Modifed on:
//'  
//'  October 01, 2016 - Created date, adapated from lda_fgs_st.cpp 
//'
// [[Rcpp::export]]
List lda_acgs_st(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    arma::mat h_grid,
    arma::mat st_grid,
    List st_grid_nbrs,
    unsigned int init_st_grid_index,
    arma::vec zetas,
    unsigned int tuning_iter,
    unsigned int max_iter_tuning,
    unsigned int max_iter_final,
    unsigned int burn_in,
    unsigned int spacing,
    double test_set_share,
    bool save_beta,
    bool save_theta,
    bool save_lp,
    bool save_hat_ratios,
    bool save_tilde_ratios,
    int verbose
) {
  
  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int valid_samples = ceil((max_iter_final - burn_in) / (double) spacing);
  unsigned int d, i, k, iter, c, word_id, word_count, topic_id, new_topic_id, g;
  unsigned int n_d; // number of words in document d
  unsigned int num_words = 0; // number of words in the corpus
  unsigned int titer; 
  
  vector < vector < unsigned int > > doc_word_ids;
  vector < vector < unsigned int > > doc_word_zids;
  vector < vector < unsigned int > > doc_word_class;
  
  arma::vec prob;
  arma::uvec doc_word_counts = arma::zeros<arma::uvec>(num_docs); // doc lengths
  arma::mat beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  arma::mat beta_t = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_t = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  arma::vec corpus_topic_counts = arma::zeros<arma::vec>(num_topics); // corpus-level topic counts
  arma::vec perplexity;
  arma::vec log_posterior;
  arma::cube theta_samples;
  arma::cube beta_samples;
  arma::uvec st_grid_index;  
  arma::vec m_hat;
  arma::mat m_hat_ratios; 
  arma::vec m_tilde;
  arma::mat m_tilde_ratios; 
  arma::vec st_grid_m_hat; 
  
  if (save_theta) 
    theta_samples = arma::cube(num_topics, num_docs, valid_samples);
  else
    theta_samples = arma::cube(num_topics, num_docs, 1);
  
  if (save_beta) 
    beta_samples = arma::cube(num_topics, vocab_size, valid_samples);
  else    
    beta_samples = arma::cube(num_topics, vocab_size, 1);

  if (save_lp) { log_posterior = arma::zeros<arma::vec>(valid_samples); }
  if (save_hat_ratios){ m_hat_ratios = zeros<mat>(h_grid.n_cols, valid_samples); }
  if (save_tilde_ratios){ m_tilde_ratios = zeros<mat>(h_grid.n_cols, valid_samples);}
  if (test_set_share > 0){perplexity = arma::zeros<arma::vec>(valid_samples);}
  
  cout << endl << endl;
  
  if (verbose > 1){
    cout << "lda_fgs_st (c++): Number of saved samples - " << valid_samples << endl;
  }

  if (verbose > 0){
    cout << "lda_fgs_st (c++): Initializes variables and count statistics....";
  }

  assert(init_st_grid_index < st_grid.n_cols);
  
  double alpha_h = st_grid(0, init_st_grid_index); 
  double eta_h = st_grid(1, init_st_grid_index); 
  
  // document-level topic mixture is used to initialize count statistics
  arma::vec alpha_vec = arma::zeros<arma::vec>(num_topics);
  alpha_vec.fill(alpha_h);
  
  assert(test_set_share < 1.);
  assert(test_set_share >= 0.);
  
  // Calculates the document word indices
  
  for (d = 0; d < num_docs; d++){
    
    arma::umat document = as<arma::umat>(docs_tf(d));
    vector < unsigned int > word_ids;
    vector < unsigned int > word_zids;
    vector < unsigned int > word_class;
    
    for (c = 0; c < document.n_cols; c++){
      word_id = document(0,c);
      word_count = document(1,c);
      
      for (i = 0; i < word_count; i++){
        // samples z for each word
        arma::vec theta_c = sample_dirichlet(num_topics, alpha_vec);
        topic_id = sample_multinomial(theta_c);
        
        word_zids.push_back(topic_id);
        word_ids.push_back(word_id);
        word_class.push_back(0); // train word
        num_words++; // increments number of words in the corpus
      }
      doc_word_counts(d) += word_count; // increments doc word counts
    }
    
    // random selection of test words
    if (test_set_share > 0){
      n_d = doc_word_counts(d); // the document length
      arma::uvec rp_d = randperm(n_d); // gets random permutations
      unsigned int num_train_words = (1. - test_set_share) * n_d;
      for (i = num_train_words; i < n_d; i++){
        word_class[rp_d(i)] = 1; // test word
      }
    }

    // doc_word_indices.push_back(word_indices);
    doc_word_ids.push_back(word_ids);
    doc_word_zids.push_back(word_zids);
    doc_word_class.push_back(word_class);
  }
  

  
  //////////////////////////////////////////////////////////////////////////////
  // updates count statistics for training words
  //////////////////////////////////////////////////////////////////////////////
  unsigned int num_test_words = 0;
  for (d = 0; d < num_docs; d++){
    
    for (i = 0; i < doc_word_counts(d); i++){
      
      topic_id = doc_word_zids[d][i];
      word_id = doc_word_ids[d][i];
      
      if (doc_word_class[d][i] == 0){ // train words
        corpus_topic_counts(topic_id) += 1.;
        theta_counts(topic_id, d) += 1.;
        beta_counts(topic_id, word_id) += 1.;
      }
      else { // test words
        num_test_words++;
      }
      
    }
    
  }
  //////////////////////////////////////////////////////////////////////////////
  
  arma::vec pred_likelihood = arma::zeros<arma::vec>(num_test_words);
  
  
  // Gets the neighbours of each point in st.grid, from the R list object 
  
  vector < vector < unsigned int > > st_grid_nbr_indices;
  for (g = 0; g < st_grid.n_cols; g++){
    vector < unsigned int > nbr_indices; 
    vec neighbors = as<vec>(st_grid_nbrs(g));
    for (i = 0; i < neighbors.n_elem; i++){
      nbr_indices.push_back(neighbors(i));
    }
    st_grid_nbr_indices.push_back(nbr_indices);
  }

  if (verbose > 0){
    cout << "DONE." << endl;
  }
  
  if (verbose > 1){
    cout << "lda_fgs_st (c++): Number of docs: " << num_docs << endl;
    cout << "lda_fgs_st (c++): Number of total words: " << num_words << endl;
    cout << "lda_fgs_st (c++): Number of test words: " << num_test_words << endl;
    cout << "lda_fgs_st (c++): Number of topics: " << num_topics << endl;
    cout << "lda_fgs_st (c++): Vocabulary size: " << vocab_size << endl;
  }
  
  
  
  // ***************************************************************************
  // BEGIN: Serial Tempering 
  // ***************************************************************************
  
  vector <unsigned int> curr_h_nbr_indices; 
  vector <unsigned int> prop_h_nbr_indices; 
  unsigned int curr_h_index = init_st_grid_index; 
  unsigned int prop_h_index; 
  unsigned int ss_idx;
  unsigned int num_gibbs_iter;
  double pd_ratio, prior_ratio, zeta_ratio, r, num_accept; 
  double M_SMOOTHING_CONST = 1e-20;
  double st_grid_nc; // stores st-grid normalizing constant for an iteration  
  double doc_denom; 
  vec h_grid_mh; // stores h-grid \hat{M}(h) ratios for an iteration 
  vec h_grid_mt; // stores h-grid \tilde{M}(h) ratios for an iteration
  vec st_grid_mh; // stores st-grid \hat{M}(h) ratios for an iteration  
  vec st_grid_occupancy;
  mat st_grid_occupancies = zeros<mat>(st_grid.n_cols, tuning_iter); 
  mat st_grid_zetas = zeros<mat>(zetas.n_elem, tuning_iter);
  vec st_grid_pr = zeros<vec>(st_grid.n_cols); // stores prior ratios for st-grid 
  vec h_grid_pr = zeros<vec>(h_grid.n_cols); // stores prior ratios for h-grid 
  int mh_inf_count;
  bool save_flag;
  arma::vec m_hat_prev; 
  arma::vec si_m_hat = zeros<arma::vec>(h_grid.n_cols);
  arma::vec m_tilde_prev; 
  arma::vec si_m_tilde = zeros<arma::vec>(h_grid.n_cols);
  unsigned int batch_size = sqrt(valid_samples); 
  unsigned int num_batches = valid_samples / batch_size; 
  unsigned int batch_idx = 0, bi = 0; 
  mat batch_m_hat_means = zeros<mat>(h_grid.n_cols, num_batches); 
  vec batch_m_hat_mean; 
  mat batch_m_tilde_means = zeros<mat>(h_grid.n_cols, num_batches); 
  vec batch_m_tilde_mean; 
  vec m_hat_mcse; 
  vec m_tilde_mcse;
  
  for (titer = 0; titer < tuning_iter; titer++){ 
    
    cout << "lda_fgs_st (c++): Serial Tempering Tuning #" << titer + 1 << endl;
    
    st_grid_zetas.col(titer) = zetas; // saves zeta for a ST iteration  
    
    // we need to reset these after each tuning iteration 
    st_grid_occupancy = zeros<vec>(st_grid.n_cols); 
    st_grid_m_hat = zeros<vec>(st_grid.n_cols); // stores \hat{M}(h)
    m_hat = zeros<vec>(h_grid.n_cols); // stores \hat{M}(h)
    m_tilde = zeros<vec>(h_grid.n_cols); // stores \tilde{M}(h)
    ss_idx = 0; 
    num_accept = 0; 
    mh_inf_count = 0; 
    
    num_gibbs_iter = (((titer+1) == tuning_iter) ? max_iter_final : max_iter_tuning);
    
    cout << "lda_fgs_st (c++): Gibbs sampling (iterations = " << num_gibbs_iter << ")" << endl;
    
    // *************************************************************************
    // BEGIN: Gibbs Sampling Loop 
    // *************************************************************************

    for (iter = 0; iter < num_gibbs_iter; iter++) { // for each Gibbs iteration
      
      if (verbose > 1) { cout << "lda_fgs_st (c++): gibbs iter #" << iter + 1; }

      save_flag = (iter >= burn_in) && (iter % spacing == 0);
      alpha_h = st_grid(0, curr_h_index); // alpha from MH 
      eta_h = st_grid(1, curr_h_index); // alpha from MH 
      
      // samples \beta
      for(k = 0; k < num_topics; k++)
        beta_t.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);

      
      for (d = 0; d < num_docs; d++) { // for each document
        
        // samples \theta
        theta_t.col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + alpha_h);
        
        n_d = doc_word_counts(d); // number of words in document d
        doc_denom = n_d - 1. + alpha_h * num_topics; // it's a constant for a term
        vector < unsigned int > word_ids = doc_word_ids[d];
        vector < unsigned int > word_zids = doc_word_zids[d];
        vector < unsigned int > word_class = doc_word_class[d];
        
        for (i = 0; i < n_d; i++) { // for each word
          
          if (word_class[i]) continue; // ignores test words from sampling
          
          topic_id = word_zids[i];
          word_id = word_ids[i];
          prob = arma::zeros <arma::vec> (num_topics); // initialize probability vector
          
          // decrements the counts by one, to ignore the current sampling word
          theta_counts(topic_id, d) -= 1.;
          beta_counts(topic_id, word_id) -= 1.;
          corpus_topic_counts(topic_id) -= 1.;
          
          // samples z's, i.e., the ***new*** topic

          for (k = 0; k < num_topics; k++){ // for each topic
            // computes p(z_{jdi} == j | \bz^{(-jdi)}, \bw, \pi)
            prob(k) = (((theta_counts(k, d) + alpha_h) / doc_denom)
                         * ((beta_counts(k, word_id) + eta_h)
                              / (corpus_topic_counts(k) + eta_h * vocab_size )));
          }
          new_topic_id = sample_multinomial(prob); // the ***new*** topic
          
          
          // increments the counts by one
          theta_counts(new_topic_id, d) += 1.;
          beta_counts(new_topic_id, word_id) += 1.;
          corpus_topic_counts(new_topic_id) += 1.;
          
          // updates newly generated topic to the database
          word_zids[i] = new_topic_id;
          
        } // for each word
        
        doc_word_zids[d] = word_zids; // updates global variable
        
      } // for each document
      
      

      if (save_flag){ // Handles burn in period
        
        // *********************************************************************
        // DO NOT CHANGE THE ORDER OF THIS CODE BLOCK 
        // Dependency: h_grid_mh calculation
        // *********************************************************************
        
        // Calculates \nu_h(\psi) / \nu_{h_(i)}(\psi) for st-grid  
        
        // Note: We compute ratios of priors instead of actual priors because  
        // sometimes priors may become too large to handle. 
        st_grid_nc = 0.0; // normalizing constant 
        for (g = 0; g < st_grid.n_cols; g++){
          prior_ratio = calc_prior_ratio(beta_t, 
                                         theta_t, 
                                         st_grid(0, g), // alpha 
                                         st_grid(1, g), // eta 
                                         alpha_h, // st alpha 
                                         eta_h); // st eta  
          st_grid_pr(g) = prior_ratio;
          st_grid_nc += (prior_ratio / ((double)st_grid.n_cols * zetas(g)));
        }

        
        // Computes the st-grid \hat{M}(h) as an online average 
        
        assert(st_grid_nc > 0);
        st_grid_mh = st_grid_pr / st_grid_nc;
        st_grid_m_hat = ((ss_idx * st_grid_m_hat + st_grid_mh) / (ss_idx + 1.)); 
        
        
        
        // *********************************************************************
        // DO NOT CHANGE THE ORDER OF THIS CODE BLOCK 
        // *********************************************************************
        
        if ((titer + 1) == tuning_iter){ // THE LAST TUNING ITERARION 
          
          // Calculates \nu_h(\psi) / \nu_{h_(i)}(\psi) for every h in h-grid 
          // given an h in st-grid  
          int inf_count = 0;
          for (g = 0; g < h_grid.n_cols; g++){
            prior_ratio = calc_prior_ratio(beta_t, 
                                           theta_t,
                                           h_grid(0, g), 
                                           h_grid(1, g), 
                                           alpha_h, 
                                           eta_h);
            h_grid_pr(g) = prior_ratio;
            if (!is_finite(prior_ratio)) inf_count++;
          }
          if (verbose > 1 && inf_count > 0) { cout << " count(NaNs): " << inf_count; }
          
          // Computes the h-grid \hat{M}(h) as an online average 
          h_grid_mh =  h_grid_pr / st_grid_nc; //st_grid_nc is *from* st-grid 
          m_hat_prev = m_hat; // saves the previous value for computing variance 
          m_hat = ((ss_idx * m_hat + h_grid_mh) / (ss_idx + 1.));  
          
          if (ss_idx == 0){
            si_m_hat = zeros<arma::vec>(h_grid.n_cols);  
          } else {
            si_m_hat += (h_grid_mh - m_hat_prev) % (h_grid_mh - m_hat);
          }
          
          // Computes the h-grid \tilde{M}(h) as an online average 
          h_grid_mt = h_grid_pr * zetas(curr_h_index); 
          m_tilde_prev = m_tilde; // saves the previous value for computing variance 
          m_tilde = ((ss_idx * m_tilde + h_grid_mt) / (ss_idx + 1.));  
          
          if (ss_idx == 0){
            si_m_tilde = zeros<arma::vec>(h_grid.n_cols);  
          } else {
            si_m_tilde += (h_grid_mt - m_tilde_prev) % (h_grid_mt - m_tilde);
          }
          
          // ///////////////////////////////////////////////////////////////////
          // Computes batch means for both tilde & hat estimates 
          // ///////////////////////////////////////////////////////////////////
          
          if (mod(ss_idx, batch_size) == 0){ // resets mean 
            if (bi > 0){ // not for the first iteration  
              batch_m_hat_means.col(batch_idx) = batch_m_hat_mean; // saves the batch mean  
              batch_m_tilde_means.col(batch_idx) = batch_m_tilde_mean; // saves the batch mean  
              batch_idx += 1; // increments batch index 
            }
            batch_m_hat_mean = h_grid_mh;  
            batch_m_tilde_mean = h_grid_mt;  
            bi = 1; // resets batch elements index 
          } else {
            batch_m_hat_mean = (bi * batch_m_hat_mean + h_grid_mh) / (bi + 1.);  
            batch_m_tilde_mean = (bi * batch_m_tilde_mean + h_grid_mt) / (bi + 1.);  
            bi += 1; // increments batch elements index 
          }

          
          // saves beta  
          if (save_beta) 
            beta_samples.slice(ss_idx) = beta_t;
          else if ((iter + 1) == max_iter_final)
            beta_samples.slice(0) = beta_t;

          // saves \theta
          if (save_theta) 
            theta_samples.slice(ss_idx) = theta_t;
          else if ((iter + 1) == max_iter_final) 
            theta_samples.slice(0) = theta_t;
          
          // saves h-grid \hat{M} ratios
          if (save_hat_ratios)
            m_hat_ratios.col(ss_idx) = h_grid_mh; 
          
          // saves h-grid \tilde{M} ratios
          if (save_tilde_ratios)
            m_tilde_ratios.col(ss_idx) = h_grid_mt;
          
          
          
          //////////////////////////////////////////////////////////////////////////
          // Computes the estimate of the predictive likelihood
          // p(w^{\text{test}}_{jdi} | \bw^{\text{train}})
          // for each test (held-out) word w^{\text{test}}_{jdi} in the corpus and 
          // the perplexity of the held-out set
          //////////////////////////////////////////////////////////////////////////
          if (test_set_share > 0){
            arma::mat beta_hat = beta_counts + eta_h; // K x V matrix
            beta_hat.each_col() /= (corpus_topic_counts + eta_h * vocab_size); // 1 x K vector
            
            unsigned int test_word_idx = 0;
            for (d = 0; d < num_docs; d++) { // for each document
              
              n_d = doc_word_counts(d); // number of words in document d
              vector < unsigned int > word_ids = doc_word_ids[d];
              vector < unsigned int > word_class = doc_word_class[d];
              
              arma::vec theta_d_hat = theta_counts.col(d) + alpha_h;
              theta_d_hat /= (n_d + alpha_h * num_topics);
              
              for (i = 0; i < n_d; i++) { // for each word
                if (word_class[i] == 0) continue; // only for test words
                
                // Calculates the predictive likelihood via online average
                double pll_t = arma::sum(theta_d_hat % beta_hat.col(word_ids[i]));
                pred_likelihood(test_word_idx) = ((ss_idx * pred_likelihood(test_word_idx) + pll_t) / (ss_idx + 1.));
                
                test_word_idx++;
              }
              
            }
            
            // perplexity of the test set, for the current iteration
            perplexity(ss_idx) = exp(-arma::mean(arma::log(pred_likelihood)));
            
            if (verbose > 1){ cout << " perp: " << perplexity(ss_idx); }
          }
          //////////////////////////////////////////////////////////////////////////
          
          if (save_lp){
            double logp = lda_log_posterior(doc_word_counts,
                                            theta_t,
                                            beta_t,
                                            doc_word_ids,
                                            doc_word_zids,
                                            alpha_h,
                                            eta_h);
            log_posterior(ss_idx) = logp;
            if (verbose > 1){ cout << " lp: " << logp; }
          }
          
        }  // THE LAST TUNING ITERARION 
        
        ss_idx += 1;
        
      } // Handles burn in period


      // ***********************************************************************
      // BEGIN: Metropolis Hastings/Serial Tempering Jump
      // ***********************************************************************
      
      st_grid_occupancy(curr_h_index) += 1.; // updates helper-grid occupancy
      
      curr_h_nbr_indices = st_grid_nbr_indices[curr_h_index];
      prop_h_index = curr_h_nbr_indices[sample_uniform_int(curr_h_nbr_indices.size())];
      prop_h_nbr_indices = st_grid_nbr_indices[prop_h_index];
      
      // Computing the acceptance ratio 
      pd_ratio = ((double)curr_h_nbr_indices.size() / (double)prop_h_nbr_indices.size()); 
      prior_ratio = calc_prior_ratio(beta_t, 
                                     theta_t, 
                                     st_grid(0, prop_h_index), // proposed alpha 
                                     st_grid(1, prop_h_index), // proposed eta
                                     alpha_h, // current alpha 
                                     eta_h); // current eta 
      
      // TODO: sometimes the prior ratio can be 0 or infinity
      // How do we handle it? Not sure yet!!!! 
      if (is_finite(prior_ratio)){
        
        zeta_ratio = zetas(curr_h_index) / zetas(prop_h_index); 
        r = min(1.0000, pd_ratio * prior_ratio * zeta_ratio);
        if (verbose > 1) {
          printf(" (r=%.3f,", r);
          printf(" pdr=%.3f,", pd_ratio); // proposal density ratio 
          printf(" ppr=%.3f,", zeta_ratio); // the ratio of normalizing constants  
          printf(" pr=%.3f)", prior_ratio); // Bayes factor ratio 
        }
        if (runif(1)(0) < r){ // acceptance 
          curr_h_index = prop_h_index; // updates the current h index  
          
          if (verbose > 1) 
            cout << " A"; 
          else if (verbose > 0)
            cout << "+"; 

          num_accept += 1.; 
        }
        else {
          if (verbose > 1) 
            cout << " R"; 
          else if (verbose > 0)
            cout << "-"; 
        }
        
      }
      else { 
        // Conditional acceptance --- This is a HACK
        // curr_h_index = prop_h_index; // updates the current h index  
        if (verbose > 1) 
          cout << " R*"; 
        else if (verbose > 0)
          cout << "*"; 
        // num_accept += 1.;       
        mh_inf_count++; 
      }
      
      // ***********************************************************************
      // END: Metropolis Hastings/Serial Tempering Jump
      // ***********************************************************************
      
      if (verbose > 1) cout << endl;

      
    } // for each Gibbs iteration
    
    // *************************************************************************
    // END: Gibbs Sampling Loop 
    // *************************************************************************
    
    if (((titer + 1) == tuning_iter) && (ss_idx >= 2)) { // THE LAST TUNING ITERARION 
      // Computes MCSE of \hat{M}(h) via Cosistant Batch Means (Flegal et al. 2008) 
      batch_m_hat_means.each_col() -= m_hat; 
      batch_m_hat_means %= batch_m_hat_means;
      m_hat_mcse = arma::sqrt(((double)batch_size / (double)(num_batches - 1)) * arma::sum(batch_m_hat_means, 1)); // the sum of elements in each row 
      m_hat_mcse /= sqrt(valid_samples); // via CBM 
      
      // Computes MCSE of \tilde{M}(h) via Cosistant Batch Means (Flegal et al. 2008) 
      batch_m_tilde_means.each_col() -= m_tilde; 
      batch_m_tilde_means %= batch_m_tilde_means; 
      m_tilde_mcse = arma::sqrt(((double)batch_size / (double)(num_batches - 1)) * arma::sum(batch_m_tilde_means, 1)); // the sum of elements in each row 
      m_tilde_mcse /= sqrt(valid_samples); // via CBM 
      
      // Computes standard deviation from running variance 
      si_m_hat = arma::sqrt(si_m_hat / (ss_idx - 1)); 
      si_m_tilde = arma::sqrt(si_m_tilde / (ss_idx - 1)); 
    }
    
    cout << endl;
    cout << "lda_fgs_st (c++): End of Gibbs sampling" << endl;
    cout << "lda_fgs_st (c++): Number of saved samples = " << ss_idx << endl;
    cout << "lda_fgs_st (c++): Serial Tempering acceptance (%) = " << 100. * (num_accept / num_gibbs_iter) << endl;
    cout << "lda_fgs_st (c++): Number of infinite prior ratios in MH Jump = " << mh_inf_count << endl;
    
    // *************************************************************************
    // BEGIN: Serial Tempering Tuning 
    // *************************************************************************
    
    st_grid_occupancies.col(titer) = st_grid_occupancy;
    cout << "lda_fgs_st (c++): Serial Tempering tuning zetas..." << endl;
    
    // Smoothing the values of \hat{M} by a small constant, M_SMOOTHING_CONST
    uvec st_grid_nf_idx = find_nonfinite(st_grid_m_hat);
    cout << "lda_fgs_st (c++): Number of nonfinite st-grid M values = " << st_grid_nf_idx.n_elem << endl;
    st_grid_m_hat(st_grid_nf_idx).zeros();
    st_grid_m_hat += M_SMOOTHING_CONST;
    
    // Note: We do tuning based on \hat{M} estimates. Another option is to tune 
    // using tilde estimates. 
    // \hat{M}(h) / \hat{M}(h1) = zeta_h_j 
    // \tilde{M}(h) / \tilde{M}(h1) = zeta_h_j 
    zetas = st_grid_m_hat / st_grid_m_hat(init_st_grid_index); 
    cout << "min(zetas) = " << zetas.min() << endl;
    cout << "max(zetas) = " << zetas.max() << endl;
    cout << "DONE." << endl << endl;
    
    // *************************************************************************
    // END: Serial Tempering Tuning 
    // *************************************************************************
    
  }
  // ***************************************************************************
  // END: Serial Tempering 
  // ***************************************************************************
  
  return List::create(
//     Named("corpus_topic_counts") = wrap(corpus_topic_counts),
//     Named("theta_counts") = wrap(theta_counts),
//     Named("beta_counts") = wrap(beta_counts),
    Named("beta_samples") = wrap(beta_samples),
    Named("theta_samples") = wrap(theta_samples),
    Named("log_posterior") = wrap(log_posterior),
    Named("perplexity") = wrap(perplexity), 
    Named("m_hat_ratios") = wrap(m_hat_ratios),
    Named("m_tilde_ratios") = wrap(m_tilde_ratios),
    Named("m_hat") = wrap(m_hat),
    Named("m_tilde") = wrap(m_tilde),
    Named("m_hat_mcse") = wrap(m_hat_mcse),
    Named("m_tilde_mcse") = wrap(m_tilde_mcse),
    Named("st_grid_zetas") = wrap(st_grid_zetas), 
    Named("st_grid_occupancies") = wrap(st_grid_occupancies), 
    Named("sigma_hat") = wrap(si_m_hat),
    Named("sigma_tilde") = wrap(si_m_tilde)
  );
  
}

