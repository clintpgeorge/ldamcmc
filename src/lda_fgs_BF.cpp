# include "utils.h"
# include "lda.h"


//' LDA: Estimate Bayes Factors using Full Gibbs Sampler
//'
//' Implements the Full Gibbs sampler for the LDA model---a Markov chain on
//' \eqn{(\beta, \theta, z)}. To compute perplexity, it first
//' partitions each document in the corpus into two sets of words: (a) a test
//' set (held-out set) and (b) a training set, given a user defined
//' \code{test_set_share}. Then, it runs the Markov chain based on the training
//' set and computes perplexity for the held-out set.
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_tf A list of corpus documents read from the Blei corpus using
//'             \code{\link{read_docs}} (term indices starts with 0)
//' @param alpha_h Hyperparameter for \eqn{\theta} sampling
//' @param eta_h Smoothing parameter for the \eqn{\beta} matrix
//' @param h_grid Grid of \eqn{(\alpha, \eta)} values
//' @param max_iter Maximum number of Gibbs iterations to be performed
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
List lda_fgs_BF_perplexity(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    double alpha_h,
    double eta_h,
    arma::mat h_grid,
    unsigned int max_iter,
    unsigned int burn_in,
    unsigned int spacing,
    bool save_theta,
    bool save_beta,
    bool save_lp,
    bool save_BF,
    int verbose,
    double test_set_share
) {

  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  unsigned int d, i, k, iter, c, word_id, word_count, topic_id, new_topic_id, g;
  unsigned int n_d; // number of words in document d
  unsigned int num_words = 0; // number of words in the corpus

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
  uvec h_grid_BF_sort_idx;
  vec h_grid_BF_mean = zeros<vec>(h_grid.n_cols);
  vec h_grid_BF_t = zeros<vec>(h_grid.n_cols);
  mat h_grid_BF;

  if (test_set_share > 0){
    perplexity = arma::zeros<arma::vec>(valid_samples);
  }

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


  if (save_BF) {
    h_grid_BF = zeros<mat>(h_grid.n_cols, valid_samples);
  }



  cout << endl << endl;
  if (verbose > 1){
    cout << "lda_fgs (c++): Number of saved samples - "
         << valid_samples << endl;
  }


  if (verbose > 0){
    cout << "lda_fgs (c++): Initializes variables and count statistics....";
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
      }
      else { // test words
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
    cout << "lda_fgs (c++): Number of docs: " << num_docs << endl;
    cout << "lda_fgs (c++): Number of total words: " << num_words << endl;
    cout << "lda_fgs (c++): Number of test words: " << num_test_words << endl;
    cout << "lda_fgs (c++): Number of topics: " << num_topics << endl;
    cout << "lda_fgs (c++): Vocabulary size: " << vocab_size << endl;
    cout << "lda_fgs (c++): alpha_h: " << alpha_h << endl;
    cout << "lda_fgs (c++): eta_h: " << eta_h << endl;
  }


  // The Gibbs sampling loop
  //

  if (verbose > 0){
    cout << "lda_fgs (c++): Full Gibbs sampling..." << endl;
  }
  unsigned int ss_idx = 0;


  for (iter = 0; iter < max_iter; iter++) { // for each Gibbs iteration

    if (verbose > 1) { cout << "lda_fgs (c++): gibbs iter# " << iter + 1; }
    else if (verbose > 0) { cout << "."; }

    bool save_flag = (iter >= burn_in) && (iter % spacing == 0);

    // samples \beta
    for(k = 0; k < num_topics; k++)
      beta_t.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);



    for (d = 0; d < num_docs; d++) { // for each document

      // samples \theta
      theta_t.col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + alpha_h);

      n_d = doc_word_counts(d); // number of words in document d
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

        // samples z's, the ***new*** topic

        new_topic_id = sample_multinomial(theta_t.col(d) % beta_t.col(word_id));


        // increments the counts by one
        theta_counts(new_topic_id, d) += 1.;
        beta_counts(new_topic_id, word_id) += 1.;
        corpus_topic_counts(new_topic_id) += 1.;

        // updates newly generated topic to the database
        word_zids[i] = new_topic_id;

      } // for each word

      doc_word_zids[d] = word_zids; // updates global variable

    } // for each document




    if (save_flag){

      // saves beta
      if (save_beta) {
        beta_samples.slice(ss_idx) = beta_t;
      }
      else if ((iter + 1) == max_iter) {
        beta_samples.slice(0) = beta_t;;
      }

      // saves \theta
      if (save_theta) {
        theta_samples.slice(ss_idx) = theta_t;
      }
      else if ((iter + 1) == max_iter) {
        theta_samples.slice(0) = theta_t;
      }

      //////////////////////////////////////////////////////////////////////////
      // Computes the estimate of the predictive likelihood
      // p(w^{\text{test}}_{jdi} | \bw^{\text{train}})
      // for each test (held-out) word w^{\text{test}}_{jdi} in the corpus and
      // the perplexity of the held-out set
      //////////////////////////////////////////////////////////////////////////
      if (test_set_share > 0) {
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

      // Computing average of Bayes Factor ratios via incremental average.
      for (g = 0; g < h_grid.n_cols; g++){
        h_grid_BF_t(g) = calc_prior_ratio(beta_t, theta_t, h_grid(0, g),
                    h_grid(1, g), alpha_h, eta_h);
      }
      if (save_BF){
        h_grid_BF.col(ss_idx) = h_grid_BF_t;
      }
      h_grid_BF_mean = (ss_idx * h_grid_BF_mean + h_grid_BF_t) / (ss_idx + 1);


      if (verbose > 1) {
        h_grid_BF_sort_idx = sort_index(h_grid_BF_mean, "descend");
        cout << " h: (" << h_grid(0, h_grid_BF_sort_idx(0)); // shows the max value
        cout << ", " << h_grid(1, h_grid_BF_sort_idx(0));
        cout << ") avg. B(h, h*):" << h_grid_BF_mean(h_grid_BF_sort_idx(0));
      }

    }

    if (verbose > 1) { cout << endl; }


    if (save_flag) {
      ss_idx += 1;
    }

  } // for each Gibbs iteration

  if (verbose > 0){
    cout << endl;
    cout << "lda_fgs (c++): Completed sampling." << endl;
  }


  return List::create(Named("corpus_topic_counts") = wrap(corpus_topic_counts),
                      Named("theta_counts") = wrap(theta_counts),
                      Named("beta_counts") = wrap(beta_counts),
                      Named("beta_samples") = wrap(beta_samples),
                      Named("theta_samples") = wrap(theta_samples),
                      Named("log_posterior") = wrap(log_posterior),
                      Named("perplexity") = wrap(perplexity),
                      Named("h_grid_BF_mean") = wrap(h_grid_BF_mean),
                      Named("h_grid_BF") = wrap(h_grid_BF)
  );

}