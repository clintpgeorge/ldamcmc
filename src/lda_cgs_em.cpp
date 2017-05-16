# include "utils.h"
# include "lda.h"


//' LDA: Gibbs-EM with Perplexity Computation
//'
//' This implements the Gibbs-EM algorithm for LDA that is mentioned in the  
//' paper Topic Modeling: Beyond Bag-of-Words. Wallach (2006). 
//' 
//' It uses the LDA collapsed Gibbs sampler---a Markov chain on \eqn{z} for the
//' E-step, and Minka (2003) fixed point iterations to optimize \eqn{h = (\eta,
//' \alpha)} in the M-step. To compute perplexity, it first partitions each
//' document in the corpus into two sets of words: (a) a test set (held-out set)
//' and (b) a training set, given a user defined \code{test_set_share}. Then, it
//' runs the Markov chain based on the training set and computes perplexity for
//' the held-out set.
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_tf A list of corpus documents read from the Blei corpus using
//'             \code{\link{read_docs}} (term indices starts with 0)
//' @param alpha_h Hyperparameter for \eqn{\theta} sampling
//' @param eta_h Smoothing parameter for the \eqn{\beta} matrix
//' @param em_max_iter Maximum number of EM iterations to be performed
//' @param gibbs_max_iter Maximum number of Gibbs iterations to be performed
//' @param burn_in Burn-in-period for the Gibbs sampler
//' @param spacing Spacing between the stored samples (to reduce correlation)
//' @param save_beta if 0 the function does not save \eqn{\beta} samples
//' @param save_theta if 0 the function does not save \eqn{\theta} samples
//' @param save_lp if 0 the function does not save computed log posterior for
//'                iterations
//' @param verbose from {0, 1, 2}
//' @param test_set_share proportion of the test words in each document. Must be
//'                       between 0. and 1.
//'
//' @return The Markov chain output as a list of
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
// [[Rcpp::export]]
List lda_cgs_em_perplexity(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    double alpha_h,
    double eta_h,
    unsigned int em_max_iter, 
    unsigned int gibbs_max_iter,
    unsigned int burn_in,
    unsigned int spacing,
    bool save_theta,
    bool save_beta,
    bool save_lp,
    int verbose,
    double test_set_share
) {
  
  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int valid_samples = ceil((gibbs_max_iter - burn_in) / (double) spacing);
  unsigned int d, i, k, iter, c, word_id, word_count, topic_id, new_topic_id;
  unsigned int n_d; // number of words in document d
  unsigned int num_words = 0; // number of words in the corpus
  double doc_denom;
  
  vector < vector < unsigned int > > doc_word_ids;
  vector < vector < unsigned int > > doc_word_zids;
  vector < vector < unsigned int > > doc_word_class;
  
  arma::uvec doc_word_counts = arma::zeros<arma::uvec>(num_docs); // doc lengths
  arma::vec corpus_topic_counts = arma::zeros<arma::vec>(num_topics); // corpus-level topic counts
  arma::vec prob;
  arma::vec perplexity; 
  arma::vec log_posterior;
  arma::cube theta_samples;
  arma::cube beta_samples;
  arma::mat beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  vector < double > opt_alphas; 
  vector < double > opt_etas; 
  
  if (save_theta) {
    theta_samples = arma::cube(num_topics, num_docs, valid_samples);
  } else {
    theta_samples = arma::cube(num_topics, num_docs, 1);
  }
  if (save_beta) {
    beta_samples = arma::cube(num_topics, vocab_size, valid_samples);
  } else {
    beta_samples = arma::cube(num_topics, vocab_size, 1);
  }
  if (save_lp) { 
    log_posterior = arma::zeros<arma::vec>(valid_samples); 
  }
  if (test_set_share > 0){
    perplexity = arma::zeros<arma::vec>(valid_samples);
  }
  
  cout << endl << endl;
  if (verbose > 1){
    cout << "lda_cgs_em (c++): Number of saved samples - " 
         << valid_samples << endl;
  }
  
  if (verbose > 0){
    cout << "lda_cgs_em (c++): Initializes variables and count statistics....";
  }
  
  // corpus level topic mixture is used to initialize count statistics
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
      } else { // test words
        num_test_words++;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  
  arma::vec pred_likelihood = arma::zeros<arma::vec>(num_test_words);
  
  if (verbose > 0){
    cout << "DONE" << endl;
  }
  
  if (verbose > 1){
    cout << "lda_cgs_em (c++): Number of docs: " << num_docs << endl;
    cout << "lda_cgs_em (c++): Number of total words: " << num_words << endl;
    cout << "lda_cgs_em (c++): Number of test words: " << num_test_words << endl;
    cout << "lda_cgs_em (c++): Number of topics: " << num_topics << endl;
    cout << "lda_cgs_em (c++): Vocabulary size: " << vocab_size << endl;
    cout << "lda_cgs_em (c++): alpha_h: " << alpha_h << endl;
    cout << "lda_cgs_em (c++): eta_h: " << eta_h << endl;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////
  // EM iterations   
  //////////////////////////////////////////////////////////////////////////////
  
  unsigned int em_iter = 0;
  double alpha_num_avg = 0; 
  double alpha_den_avg = 0; 
  double eta_num_avg = 0; 
  double eta_den_avg = 0; 
  double eta_num_sum = 0; 
  double eta_den_sum = 0; 
  double alpha_num_sum = 0; 
  double alpha_den_sum = 0; 
  double alpha_h_old; 
  double eta_h_old; 
  cout.precision(10);
  
  while ((em_iter <= 2) || (em_iter <= em_max_iter)) { // for EM iteration 

    em_iter++;
    
    ////////////////////////////////////////////////////////////////////////////
    // E Step: Gibbs sampling    
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose > 0){
      cout << "lda_cgs_em (c++): em-iter #" << em_iter << endl << endl;
      cout << "lda_cgs_em (c++): Collapsed Gibbs sampling..." << endl;
    }
    
    unsigned int ss_idx = 0;

    
    for (iter = 0; iter < gibbs_max_iter; iter++) { // for each Gibbs iteration
      
      if (verbose > 1) { cout << "lda_cgs_em (c++): gibbs iter# " << iter + 1; }
      else if (verbose > 0) { cout << "."; }
      
      bool save_flag = (iter >= burn_in) && (iter % spacing == 0);
      
      // samples \beta
      if (save_flag && save_beta) {
        for(k = 0; k < num_topics; k++)
          beta_samples.slice(ss_idx).row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);
      }
      else if (save_flag && ((iter + 1) == gibbs_max_iter)) {
        for(k = 0; k < num_topics; k++)
          beta_samples.slice(0).row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);
      }
      
      ////////////////////////////////////////////////////////////////////////////
      // updates z's
      ////////////////////////////////////////////////////////////////////////////
      
      alpha_num_sum = 0; 
      alpha_den_sum = 0; 
      
      for (d = 0; d < num_docs; d++) { // for each document
        
        // samples \theta
        if (save_flag && save_theta) {
          theta_samples.slice(ss_idx).col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + alpha_h);
        }
        else if (save_flag && ((iter + 1) == gibbs_max_iter)) {
          theta_samples.slice(0).col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + alpha_h);
        }
        
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
          
          // samples z's
          
          for (k = 0; k < num_topics; k++){ // for each topic
            // computes p(z_{jdi} == j | \bz^{(-jdi)}, \bw, \pi)
            prob(k) = (((theta_counts(k, d) + alpha_h) / doc_denom)
                         * ((beta_counts(k, word_id) + eta_h)
                              / (corpus_topic_counts(k) + eta_h * vocab_size )));
          }
          new_topic_id = sample_multinomial(prob); // the ***new*** topic
          
          // assert(new_topic_id < num_topics); // check
          
          // increments the counts by one
          theta_counts(new_topic_id, d) += 1.;
          beta_counts(new_topic_id, word_id) += 1.;
          corpus_topic_counts(new_topic_id) += 1.;
          
          // updates newly generated topic to the database
          word_zids[i] = new_topic_id;
          
        } // for each word
        
        doc_word_zids[d] = word_zids; // updates global variable
        
        // Computes \alpha statistics 
        
        alpha_num_sum += sum(digamma_vec(theta_counts.col(d) + alpha_h) - Rf_digamma(alpha_h)); 
        alpha_den_sum += num_topics * (Rf_digamma(doc_word_counts(d) + num_topics * alpha_h) - Rf_digamma(num_topics * alpha_h)); 
        
      } // for each document
      
      //////////////////////////////////////////////////////////////////////////
      // Computes \alpha and \eta statistics 
      //////////////////////////////////////////////////////////////////////////
      if (save_flag){
        alpha_num_avg = ((ss_idx * alpha_num_avg + (alpha_num_sum / (double) num_docs)) / (ss_idx + 1.)); 
        alpha_den_avg = ((ss_idx * alpha_den_avg + (alpha_den_sum / (double) num_docs)) / (ss_idx + 1.)); 
        
        eta_num_sum = 0; 
        eta_den_sum = 0; 
        for(k = 0; k < num_topics; k++){
          eta_num_sum += sum(digamma_rowvec(beta_counts.row(k) + eta_h) - Rf_digamma(eta_h)); 
          eta_den_sum += vocab_size * (Rf_digamma(corpus_topic_counts(k) + vocab_size * eta_h) - Rf_digamma(vocab_size * eta_h)); 
        }
        eta_num_avg = ((ss_idx * eta_num_avg + (eta_num_sum / (double) num_topics)) / (ss_idx + 1.)); 
        eta_den_avg = ((ss_idx * eta_den_avg + (eta_den_sum / (double) num_topics)) / (ss_idx + 1.)); 
      }
      
      
      
      ////////////////////////////////////////////////////////////////////////////
      // Computes the estimate of the predictive likelihood
      // p(w^{\text{test}}_{jdi} | \bw^{\text{train}})
      // for each test (held-out) word w^{\text{test}}_{jdi} in the corpus and the
      // perplexity of the held-out set
      ////////////////////////////////////////////////////////////////////////////
      if (save_flag && test_set_share > 0){
        
        arma::mat beta_hat = beta_counts + eta_h; // K x V matrix
        beta_hat.each_col() /= (corpus_topic_counts + eta_h * vocab_size); // arma::sum(beta_hat, 1); // 1 x K vector
        
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
      
      
      
      if (save_flag && save_beta && save_theta && save_lp){
        double logp = lda_log_posterior(doc_word_counts,
                                        theta_samples.slice(ss_idx),
                                        beta_samples.slice(ss_idx),
                                        doc_word_ids,
                                        doc_word_zids,
                                        alpha_h,
                                        eta_h);
        log_posterior(ss_idx) = logp;
        if (verbose > 1){ cout << " lp: " << logp; }
      }
      
      
      if (verbose > 1) { cout << endl; }
      if (save_flag) { ss_idx += 1; }
      
    } // for each Gibbs iteration
    
    if (verbose > 0){
      cout << "lda_cgs_em (c++): Completed sampling." << endl;
      cout << endl;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // M Step: Hyperparameter Optimization 
    ////////////////////////////////////////////////////////////////////////////
    
    // new alpha 
    alpha_h_old = alpha_h; 
    alpha_h = alpha_h * alpha_num_avg / alpha_den_avg; 
    
    // new eta 
    eta_h_old = eta_h; 
    eta_h = eta_h * eta_num_avg / eta_den_avg; 
    
    opt_alphas.push_back(alpha_h); 
    opt_etas.push_back(eta_h); 
    
    cout << "lda_cgs_em (c++): new alpha: " << alpha_h << " old-alpha: " << alpha_h_old; 
    cout << " alpha_num_avg: " << alpha_num_avg << " alpha_den_avg: " << alpha_den_avg << endl; 
    cout << "lda_cgs_em (c++): new eta: " << eta_h << " old-eta: " << eta_h_old; 
    cout << " eta_num_avg: " << eta_num_avg << " eta_den_avg: " << eta_den_avg << endl; 
    cout << endl; 
    
  } // for EM iteration 
  
  
  
  return List::create(Named("corpus_topic_counts") = wrap(corpus_topic_counts),
                      Named("theta_counts") = wrap(theta_counts),
                      Named("beta_counts") = wrap(beta_counts),
                      Named("beta_samples") = wrap(beta_samples),
                      Named("theta_samples") = wrap(theta_samples),
                      Named("log_posterior") = wrap(log_posterior),
                      Named("perplexity") = wrap(perplexity), 
                      Named("opt_alphas") = wrap(opt_alphas),
                      Named("opt_etas") = wrap(opt_etas));
  
}

//' LDA: Gibbs-EM
//'
//' This implements the Gibbs-EM algorithm for LDA that is mentioned in the  
//' paper Topic Modeling: Beyond Bag-of-Words. Wallach (2006). 
//'
//' @export
//'
//' @family MCMC
//'
// [[Rcpp::export]]
List lda_cgs_em(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    double alpha_h,
    double eta_h,
    unsigned long int em_max_iter, 
    unsigned long int gibbs_max_iter,
    unsigned long int burn_in,
    unsigned long int spacing,
    int verbose
) {
  
  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int d, i, k, c, word_id, word_count, topic_id, new_topic_id;
  unsigned long int iter, em_iter;  
  unsigned int n_d; // number of words in document d
  unsigned int num_words = 0; // number of words in the corpus
  double doc_denom;
  
  vector < vector < unsigned int > > doc_word_ids;
  vector < vector < unsigned int > > doc_word_zids;
  
  arma::uvec doc_word_counts = arma::zeros<arma::uvec>(num_docs); // doc lengths
  arma::vec corpus_topic_counts = arma::zeros<arma::vec>(num_topics); // corpus-level topic counts
  arma::vec prob;
  arma::mat beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  vector < double > opt_alphas; 
  vector < double > opt_etas; 
  
  if (verbose > 0){
    cout << "lda_cgs_em (c++): Initializes variables and count statistics....";
  }
  
  // corpus level topic mixture is used to initialize count statistics
  arma::vec alpha_vec = arma::zeros<arma::vec>(num_topics);
  alpha_vec.fill(alpha_h);

  // Calculates the document word indices
  
  for (d = 0; d < num_docs; d++){
    
    arma::umat document = as<arma::umat>(docs_tf(d));
    vector < unsigned int > word_ids;
    vector < unsigned int > word_zids;
    
    for (c = 0; c < document.n_cols; c++){
      word_id = document(0,c);
      word_count = document(1,c);
      
      for (i = 0; i < word_count; i++){
        // samples z for each word
        arma::vec theta_c = sample_dirichlet(num_topics, alpha_vec);
        topic_id = sample_multinomial(theta_c);
        
        word_zids.push_back(topic_id);
        word_ids.push_back(word_id);
        num_words++; // increments number of words in the corpus
        
        // updates count statistics
        corpus_topic_counts(topic_id) += 1.;
        theta_counts(topic_id, d) += 1.;
        beta_counts(topic_id, word_id) += 1.;
      }
      doc_word_counts(d) += word_count; // increments doc word counts
    }
    
    // doc_word_indices.push_back(word_indices);
    doc_word_ids.push_back(word_ids);
    doc_word_zids.push_back(word_zids);
  
  }
  
  if (verbose > 0){
    cout << "DONE" << endl;
  }
  
  if (verbose > 1){
    cout << "lda_cgs_em (c++): Number of docs: " << num_docs << endl;
    cout << "lda_cgs_em (c++): Number of total words: " << num_words << endl;
    cout << "lda_cgs_em (c++): Number of topics: " << num_topics << endl;
    cout << "lda_cgs_em (c++): Vocabulary size: " << vocab_size << endl;
    cout << "lda_cgs_em (c++): Initial alpha: " << alpha_h << endl;
    cout << "lda_cgs_em (c++): Initial eta: " << eta_h << endl;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////
  // EM iterations   
  //////////////////////////////////////////////////////////////////////////////
  
  double alpha_num_avg = 0; 
  double alpha_den_avg = 0; 
  double eta_num_avg = 0; 
  double eta_den_avg = 0; 
  double eta_num_sum = 0; 
  double eta_den_sum = 0; 
  double alpha_num_sum = 0; 
  double alpha_den_sum = 0; 
  double alpha_h_old; 
  double eta_h_old; 
  cout.precision(10);
  em_iter = 0;
  
  while ((em_iter <= 2) || (em_iter <= em_max_iter)) { // for EM iteration 

    em_iter++;
    
    ////////////////////////////////////////////////////////////////////////////
    // E Step: Gibbs sampling    
    ////////////////////////////////////////////////////////////////////////////
    
    // Increments the maximum number of iterations for the Gibbs sampling run 
    // in each EM iteration (We did not find any improvement, hence commented 
    // the code.)
    // NOTE: Due to the limits of C++ integer types, 
    // see http://www.cplusplus.com/reference/climits/
    // the maximum value we can accomodate is 4294967295 (232-1)
    // unsigned long int num_iter = pow(2, (6 + em_iter));
    // unsigned long int gibbs_max = num_iter + burn_in; 
    // if (num_iter > gibbs_max_iter) {
    //   gibbs_max = gibbs_max_iter + burn_in;
    // }
    
    // a fixed maximum number of iterations
    unsigned long int gibbs_max = gibbs_max_iter; 

    if (verbose > 0){
      cout << "lda_cgs_em (c++): em-iter #" << em_iter << endl;
      // cout << "lda_cgs_em (c++): gibbs-max-iter #" << gibbs_max << " " << num_iter << endl;
      cout << "lda_cgs_em (c++): Collapsed Gibbs sampling..." << endl;
    }
    
    unsigned int ss_idx = 0;
    
    
    for (iter = 0; iter < gibbs_max; iter++) { // for each Gibbs iteration
      
      if (verbose > 1) { cout << "lda_cgs_em (c++): gibbs iter# " << iter + 1 << endl; }
      else if (verbose > 0) { cout << "."; }
      
      bool save_flag = (iter >= burn_in) && (iter % spacing == 0);
      

      ////////////////////////////////////////////////////////////////////////////
      // updates z's
      ////////////////////////////////////////////////////////////////////////////
      
      alpha_num_sum = 0; 
      alpha_den_sum = 0; 
      
      for (d = 0; d < num_docs; d++) { // for each document
        
        n_d = doc_word_counts(d); // number of words in document d
        doc_denom = n_d - 1. + alpha_h * num_topics; // it's a constant for a term
        vector < unsigned int > word_ids = doc_word_ids[d];
        vector < unsigned int > word_zids = doc_word_zids[d];

        for (i = 0; i < n_d; i++) { // for each word

          topic_id = word_zids[i];
          word_id = word_ids[i];
          prob = arma::zeros <arma::vec> (num_topics); // initialize probability vector
          
          // decrements the counts by one, to ignore the current sampling word
          theta_counts(topic_id, d) -= 1.;
          beta_counts(topic_id, word_id) -= 1.;
          corpus_topic_counts(topic_id) -= 1.;
          
          // samples z's
          
          for (k = 0; k < num_topics; k++){ // for each topic
            // computes p(z_{jdi} == j | \bz^{(-jdi)}, \bw, \pi)
            prob(k) = (((theta_counts(k, d) + alpha_h) / doc_denom)
                         * ((beta_counts(k, word_id) + eta_h)
                              / (corpus_topic_counts(k) + eta_h * vocab_size )));
          }
          new_topic_id = sample_multinomial(prob); // the ***new*** topic
          
          // assert(new_topic_id < num_topics); // check
          
          // increments the counts by one
          theta_counts(new_topic_id, d) += 1.;
          beta_counts(new_topic_id, word_id) += 1.;
          corpus_topic_counts(new_topic_id) += 1.;
          
          // updates newly generated topic to the database
          word_zids[i] = new_topic_id;
          
        } // for each word
        
        doc_word_zids[d] = word_zids; // updates global variable
        
        // Computes \alpha statistics 
        if (save_flag){
          alpha_num_sum += sum(digamma_vec(theta_counts.col(d) + alpha_h) - Rf_digamma(alpha_h)); 
          alpha_den_sum += num_topics * (Rf_digamma(doc_word_counts(d) + num_topics * alpha_h) - Rf_digamma(num_topics * alpha_h)); 
        }
        
      } // for each document
      
      //////////////////////////////////////////////////////////////////////////
      // Computes \alpha and \eta statistics as an incremental average 
      //////////////////////////////////////////////////////////////////////////
      if (save_flag){
        alpha_num_avg = ((ss_idx * alpha_num_avg + (alpha_num_sum / (double) num_docs)) / (ss_idx + 1.)); 
        alpha_den_avg = ((ss_idx * alpha_den_avg + (alpha_den_sum / (double) num_docs)) / (ss_idx + 1.)); 
        
        eta_num_sum = 0; 
        eta_den_sum = 0; 
        for(k = 0; k < num_topics; k++){
          eta_num_sum += sum(digamma_rowvec(beta_counts.row(k) + eta_h) - Rf_digamma(eta_h)); 
          eta_den_sum += vocab_size * (Rf_digamma(corpus_topic_counts(k) + vocab_size * eta_h) - Rf_digamma(vocab_size * eta_h)); 
        }
        eta_num_avg = ((ss_idx * eta_num_avg + (eta_num_sum / (double) num_topics)) / (ss_idx + 1.)); 
        eta_den_avg = ((ss_idx * eta_den_avg + (eta_den_sum / (double) num_topics)) / (ss_idx + 1.));
        
        ss_idx += 1; // increments the count 
      }

    } // for each Gibbs iteration
    
    if (verbose > 0){
      cout << "Number of saved samples: " << ss_idx << endl;
      cout << "lda_cgs_em (c++): End of Gibbs sampling." << endl;
      cout << endl;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // M Step: Hyperparameter Optimization 
    // 
    // References: 
    // 	1. Topic modeling: beyond bag of words (Wallach 2006), Eq. 30 
    // 	2. Estimating a Dirichlet distribution (Minka 2012), Eq. 55   
    ////////////////////////////////////////////////////////////////////////////
    
    // new alpha 
    alpha_h_old = alpha_h; 
    alpha_h = alpha_h * alpha_num_avg / alpha_den_avg; 
    
    // new eta 
    eta_h_old = eta_h; 
    eta_h = eta_h * eta_num_avg / eta_den_avg; 
    
    opt_alphas.push_back(alpha_h); 
    opt_etas.push_back(eta_h); 
    
    cout << "lda_cgs_em (c++): new-alpha: " << alpha_h << " old-alpha: " << alpha_h_old; 
    cout << " num (avg) / denom (avg): " << alpha_num_avg / alpha_den_avg << endl; 
    cout << "lda_cgs_em (c++): new-eta: " << eta_h << " old-eta: " << eta_h_old; 
    cout << " num (avg) / denom (avg): " << eta_num_avg / eta_den_avg << endl; 
    cout << endl; 
    
  } // end of EM iteration loop
  

  return List::create(
    Named("theta_counts") = wrap(theta_counts),
    Named("beta_counts") = wrap(beta_counts),
    Named("opt_alphas") = wrap(opt_alphas),
    Named("opt_etas") = wrap(opt_etas)
  );
  
}