# include "utils.h"
# include "lda.h"


//' LDA: Full Gibbs Sampler with Posterior Predictive Value
//'
//' Implements the Full Gibbs sampler for the LDA model---a Markov chain on
//' \eqn{(\beta, \theta, z)}. The log posterior predictive value is based on 
//' Zhe Chen (2015)
//'
//' @param num_topics Number of topics in the corpus
//' @param vocab_size  Vocabulary size
//' @param docs_tf A list of corpus documents read from the Blei corpus using
//'             \code{\link{read_docs}} (term indices starts with 0)
//' @param alpha_h Hyperparameter for \eqn{\theta} sampling
//' @param eta_h Smoothing parameter for the \eqn{\beta} matrix
//' @param max_iter Maximum number of Gibbs iterations to be performed
//' @param burn_in Burn-in-period for the Gibbs sampler
//' @param spacing Spacing between the stored samples (to reduce correlation)
//' @param verbose from {0, 1, 2}
//'
//' @return The Markov chain output as a list of
//'   \item{lppv}{log posterior predictive values of each document}
//'   \item{lppc}{averge of log posterior predictive values}
//'
//' @export
//'
//' @family MCMC
//'
// [[Rcpp::export]]
List lda_fgs_ppc(
    unsigned int num_topics,
    unsigned int vocab_size,
    List docs_tf,
    double alpha_h,
    double eta_h,
    unsigned int max_iter,
    unsigned int burn_in,
    unsigned int spacing,
    int verbose
) {

  unsigned int num_docs = docs_tf.size(); // number of documents
  unsigned int saved_samples = ceil((max_iter - burn_in) / (double) spacing);
  unsigned int d, i, k, iter, c, word_id, word_count, topic_id, new_topic_id, hod;
  unsigned int n_d; // number of words in document d
  unsigned int num_words = 0; // number of words in the corpus
  unsigned int ss_idx;

  vector < vector < unsigned int > > doc_word_ids;
  vector < vector < unsigned int > > doc_word_zids;
  vector < vector < unsigned int > > doc_word_zids_init;

  arma::vec prob;
  arma::uvec doc_word_counts = arma::zeros<arma::uvec>(num_docs); // doc lengths
  arma::mat beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  arma::mat beta_t = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
  arma::mat theta_t = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
  mat lppv = zeros<mat>(num_docs, saved_samples);
  double llw_hod, shift_hod, avg_ppv, lppv_sum = 0.0;

  cout << endl << endl;
  if (verbose > 0){
    cout << "lda_fgs (c++): Number of saved samples - "
         << saved_samples << endl;
  }
  if (verbose > 0){
    cout << "lda_fgs (c++): Initializes variables and count statistics....";
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
      }
      doc_word_counts(d) += word_count; // increments doc word counts
    }

    // doc_word_indices.push_back(word_indices);
    doc_word_ids.push_back(word_ids);
    doc_word_zids.push_back(word_zids);
    doc_word_zids_init.push_back(word_zids);
  }

  if (verbose > 0){
    cout << "DONE" << endl;
  }

  if (verbose > 1){
    cout << "lda_fgs (c++): Number of docs: " << num_docs << endl;
    cout << "lda_fgs (c++): Number of total words: " << num_words << endl;
    cout << "lda_fgs (c++): Number of topics: " << num_topics << endl;
    cout << "lda_fgs (c++): Vocabulary size: " << vocab_size << endl;
    cout << "lda_fgs (c++): alpha_h: " << alpha_h << endl;
    cout << "lda_fgs (c++): eta_h: " << eta_h << endl;
  }


  //
  // Leave One Out Posterior Predictive Value Computation:
  //

  for (hod = 0; hod < num_docs; hod++){ // For each held-out document hod

    if (verbose > 0){
      cout << "lda_fgs (c++): Gibbs sampling [doc #" << hod + 1 << "]... " << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Initializes count statistics and other variables
    ////////////////////////////////////////////////////////////////////////////
    ss_idx = 0;
    theta_counts = arma::zeros<arma::mat>(num_topics, num_docs); // K x D matrix
    beta_counts = arma::zeros<arma::mat>(num_topics, vocab_size); // K x V matrix
    for (d = 0; d < num_docs; d++){
      if (d == hod) continue; // ignores the held-out document
      for (i = 0; i < doc_word_counts(d); i++){
        topic_id = doc_word_zids[d][i] = doc_word_zids_init[d][i]; // INITIALIZES w/ the init
        word_id = doc_word_ids[d][i];
        theta_counts(topic_id, d) += 1.;
        beta_counts(topic_id, word_id) += 1.;
      }
    }
    ////////////////////////////////////////////////////////////////////////////

    for (iter = 0; iter < max_iter; iter++) { // for each Gibbs iteration

      if (verbose > 1) { cout << "\nlda_fgs (c++): gibbs iter# " << iter + 1; }
      else if (verbose > 0) { cout << "."; }

      bool save_flag = (iter >= burn_in) && (iter % spacing == 0);

      // samples \beta
      for(k = 0; k < num_topics; k++)
        beta_t.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta_h);

      for (d = 0; d < num_docs; d++) { // for each document

        if (d == hod) continue; // ignores the held-out document

        // samples \theta
        theta_t.col(d) = sample_dirichlet(num_topics, theta_counts.col(d) + alpha_h);

        n_d = doc_word_counts(d); // number of words in document d
        vector < unsigned int > word_ids = doc_word_ids[d];
        vector < unsigned int > word_zids = doc_word_zids[d];

        for (i = 0; i < n_d; i++) { // for each word

          topic_id = word_zids[i];
          word_id = word_ids[i];
          prob = arma::zeros <arma::vec> (num_topics); // initialize probability vector

          // decrements the counts by one, to ignore the current sampling word
          theta_counts(topic_id, d) -= 1.;
          beta_counts(topic_id, word_id) -= 1.;

          // samples z's, the ***new*** topic

          new_topic_id = sample_multinomial(theta_t.col(d) % beta_t.col(word_id));


          // increments the counts by one
          theta_counts(new_topic_id, d) += 1.;
          beta_counts(new_topic_id, word_id) += 1.;

          // updates newly generated topic to the database
          word_zids[i] = new_topic_id;

        } // for each word

        doc_word_zids[d] = word_zids; // updates global variable

      } // for each document

      if (save_flag) {
        // computes the log posterior predictive value for the held-out document
        // using the current iteration
        // samples \theta for held-out document hod
        // 

        vec theta_hod = sample_dirichlet(num_topics, alpha_vec);

        vector < unsigned int > word_ids = doc_word_ids[hod];
        llw_hod = 0.;
        for (i = 0; i < doc_word_counts(hod); i++){
          llw_hod += log(accu(theta_hod % beta_t.col(word_ids[i])));
        }
        lppv(hod, ss_idx) = llw_hod;
        if (verbose > 1) {
          cout << " log(PPV) = " << llw_hod;
        }

        ss_idx += 1;
      }

      if (verbose > 1) { cout << endl; }

    } // for each Gibbs iteration

    if (verbose > 0){
      cout << endl;
      cout << "lda_fgs (c++): Completed sampling." << endl;
    }

    // Computes the estimate of the predictive likelihood for
    // the held-out document 'hod'

    shift_hod = max(lppv.row(hod)) - 20;
    avg_ppv = log(mean(exp(lppv.row(hod) - shift_hod))) + shift_hod;
    lppv_sum += avg_ppv;

    if (verbose > 0){
      cout << "lda_fgs (c++): log (PPV[avg](" << hod + 1 << ")) = "
           << avg_ppv << endl;
    }

  } // For each held-out document hod


  return List::create(
    Named("lppv") = wrap(lppv),
    Named("lppc") = wrap(lppv_sum / (double)num_docs)
  );

}