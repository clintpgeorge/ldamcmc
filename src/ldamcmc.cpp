#include "ldamcmc.h"


/**
 *  Computes Log Posterior Predictive Value for the LDA model: 
 *  
 *  Implementation of the log posterior predictive value is based on Zhe Chen 
 *  (2015) [dissertation]. 
 *
 * 	Arguments:
 * 		num_topics_              - number of topics
 * 		vocab_size_              - vocabulary size
 * 		word_ids_                - vocabulary ids of each word in each corpus document
 * 		doc_lengths_             - number of words in each document as a vector
 * 		topic_assignments_       - initial topic assignment of each word in each corpus document
 * 		alpha_v_                 - hyperparameter vector for the document topic Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - max number of Gibbs iterations to be perfomed
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are stored
 *
 * 	Returns:
 * 		Log posterior predictive value of the corpus 
 *
 */
RcppExport SEXP lda_fgs_lppv(SEXP num_topics_, 
                             SEXP vocab_size_, 
                             SEXP word_ids_, 
                             SEXP doc_lengths_, 
                             SEXP topic_assignments_, 
                             SEXP alpha_v_, 
                             SEXP eta_, 
                             SEXP max_iter_, 
                             SEXP burn_in_, 
                             SEXP spacing_) {
  
  cout << "lda_fgs (c++): initialization of variables ..." << endl;
  
  // Variable from the R interface  
  
  uvec doc_lengths = as<uvec>(doc_lengths_);
  uvec word_ids = as<uvec>(word_ids_);
  uvec z = as<uvec>(topic_assignments_);
  vec alpha_v = as<vec>(alpha_v_);

  double eta = as<double>(eta_);
  unsigned int num_topics = as<unsigned int>(num_topics_);
  unsigned int vocab_size = as<unsigned int>(vocab_size_);
  unsigned int max_iter = as<unsigned int>(max_iter_);
  unsigned int burn_in = as<unsigned int>(burn_in_);
  unsigned int spacing = as<unsigned int>(spacing_);

  // Function variables 
  
  unsigned int saved_samples = ceil((max_iter - burn_in) / (double) spacing);
  unsigned int num_docs = doc_lengths.n_elem;

  mat beta_samples = zeros<mat>(num_topics, vocab_size);
  mat beta_counts = zeros<mat>(num_topics, vocab_size);
  mat lppv = zeros<mat>(num_docs, saved_samples);
  double cp, llw_hod, shift_hod, avg_ppv, lppv_sum = 0.0; 

  vector < vector < unsigned int > > doc_word_indices;
  unsigned int d, i, k, iter, ss_idx, instances = 0, hod; 
  unsigned int msg_interval = (100 >= max_iter)?1:100;
  

  // Calculates the indices for each word in a document as a vector
  
  for (d = 0; d < num_docs; d++){
    vector < unsigned int > word_idx;
    for (i = 0; i < doc_lengths(d); i++){
      word_idx.push_back(instances);
      instances++;
    }
    doc_word_indices.push_back(word_idx);
  }
  
  cout << "lda_fgs (c++): total number of words - " << instances << endl;
  cout << "lda_fgs (c++): completed initialization." << endl;  
  
  
  
  for (hod = 0; hod < num_docs; hod++){ // For each held-out document hod 

    cout << "lda_fgs (c++): Gibbs sampling [doc #" << hod+1 << "]... " << endl;
    
    ss_idx = 0; 
    
    // Initilizes beta
    
    beta_counts.fill(eta); // initializes with the smoothing parameter
    for (d = 0; d < num_docs; d++){ // for each document
      if (d == hod) continue; // ignores the held-out document    
      vector < unsigned int > word_idx = doc_word_indices[d];
      for (i = 0; i < doc_lengths(d); i++)
        beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1.0;
    }

    // Full Gibbs sampler 
    
    for (iter = 0; iter < max_iter; iter++){ // The Gibbs sampling loop
      
      if (iter % msg_interval == 0) { 
        cout << "lda_fgs (c++): gibbs iter# " << iter + 1;
      }
      
      // samples \beta
      
      for(k = 0; k < num_topics; k++)
        beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));

      for (d = 0; d < num_docs; d++){ // for each document
        
        if (d == hod) continue; // ignores the held-out document    
        
        vector < unsigned int > word_idx = doc_word_indices[d];
        
        // samples \theta
        
        vec partition_counts = alpha_v; // initializes with the smoothing parameter
        for (i = 0; i < doc_lengths(d); i++)
          partition_counts(z(word_idx[i])) += 1.0;
        vec theta_d = sample_dirichlet(num_topics, partition_counts);

        // samples z and updates \beta counts
        
        for(i = 0; i < doc_lengths(d); i++)
          beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1.0; // excludes document d's word-topic counts
        for (i = 0; i < doc_lengths(d); i++)
          z(word_idx[i]) = sample_multinomial(theta_d % beta_samples.col(word_ids(word_idx[i])));
        for(i = 0; i < doc_lengths(d); i++)
          beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1.0; // includes document d's word-topic counts
        
      }
      
      if ((iter >= burn_in) && (iter % spacing == 0)){ // the burn-in period
        
        // samples \theta for held-out document hod  
        
        vec theta_hod = sample_dirichlet(num_topics, alpha_v); 
        
        // computes the log posterior predictive value of the held-out document 
        // using the current sample/iteration  
        
        vector < unsigned int > word_idx = doc_word_indices[hod];
        llw_hod = 0.0;  
        for (i = 0; i < doc_lengths(hod); i++){
          cp = accu(theta_hod % beta_samples.col(word_ids(word_idx[i])));
          llw_hod += log(cp);
        }
        lppv(hod, ss_idx) = llw_hod;
        if (iter % msg_interval == 0) { cout << " log(PPV) = " << llw_hod; }

        ss_idx++; // updates the counter  
        
      } // the burn-in period
      
      
      if (iter % msg_interval == 0) { cout << endl; }
      
    } // The Gibbs sampling loop
    
    cout << "lda_fgs (c++): completed Gibbs sampling." << endl;
    cout << "lda_fgs (c++): number of saved samples - " << ss_idx << endl;
    
    // Computes the estimate of the predictive likelihood of the held-out 
    // document hod  
    
    shift_hod = max(lppv.row(hod)) - 20; 
    avg_ppv = log(mean(exp(lppv.row(hod) - shift_hod))) + shift_hod;
    
    cout << "lda_fgs (c++): log (PPV[avg](" << hod+1 << ")) = " << avg_ppv << endl;
    
    lppv_sum += avg_ppv;

  } // For each held-out document hod 
  
  
  return List::create(
    Named("lppv") = wrap(lppv), 
    Named("lppc") = wrap(lppv_sum/(double)num_docs)
  );
  
}



/**
 *  The LDA full Gibbs sampler (FGS): 
 *
 * 	Arguments:
 * 		num_topics_              - number of topics
 * 		vocab_size_              - vocabulary size
 * 		word_ids_                - vocabulary ids of each word in each corpus document
 * 		doc_lengths_             - number of words in each document as a vector
 * 		topic_assignments_       - initial topic assignment of each word in each corpus document
 * 		alpha_v_                 - hyperparameter vector for the document topic Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - max number of Gibbs iterations to be perfomed
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are stored
 * 		save_z_                  - save the sampled z (N x 1 vector) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_beta_         	     - save the sampled beta (K x V matrix) for each  
 * 		                           iteration, values: {1, 0} 
 * 		save_theta_              - save the sampled theta (K x D matrix) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_lp_                 - compute and save the log posterior of the LDA 
 * 		                           model
 *
 * 	Returns:
 * 		thetas                   - sampled thetas after the burn in period
 * 		betas                    - sampled betas after the burn in period
 * 		Z                        - sampled word topic assignments after the burn in period
 * 		lp                       - log posterior of the LDA model after ignoring the normalizing constants
 *
 */
RcppExport SEXP lda_fgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
                        SEXP doc_lengths_, SEXP topic_assignments_, 
                        SEXP alpha_v_, SEXP eta_, SEXP max_iter_, SEXP burn_in_, 
                        SEXP spacing_, SEXP save_z_, SEXP save_beta_, 
                        SEXP save_theta_, SEXP save_lp_) {

	cout << "lda_fgs (c++): init process..." << endl;

	// Variable from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_);
	uvec word_ids = as<uvec>(word_ids_);
	uvec z = as<uvec>(topic_assignments_);
	vec alpha_v = as<vec>(alpha_v_);

	double eta = as<double>(eta_);
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_lp = as<unsigned int>(save_lp_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = word_ids.n_elem;
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_fgs (c++): Total number of words - " << num_word_instances << endl;

	cube thetas;
	cube betas;
	umat Z; 
	vec log_posterior;
	
	if (save_z) { Z = zeros<umat>(num_word_instances, valid_samples); }
	if (save_beta) { betas = cube(num_topics, vocab_size, valid_samples); }
	if (save_theta) { thetas = cube(num_topics, num_docs, valid_samples); }
	if (save_lp) { log_posterior = zeros<vec>(valid_samples); } 

	mat beta_samples = zeros<mat>(num_topics, vocab_size);
	mat theta_samples = zeros<mat>(num_topics, num_docs);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);

	vector < vector < unsigned int > > doc_word_indices;
	unsigned int d, i, k, iter, ss_idx = 0, instances = 0;
	rowvec eta_v = zeros<rowvec>(vocab_size);
	for(k = 0; k < vocab_size; k++) { eta_v(k) = eta; }

	// Calculates the indices for each word in a document as a vector
	for (d = 0; d < num_docs; d++){
		vector < unsigned int > word_idx;
		for (i = 0; i < doc_lengths(d); i++){
			word_idx.push_back(instances);
			instances++;
		}
		doc_word_indices.push_back(word_idx);
	}

	// Initilizes beta

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1;
	}

	cout << "lda_fgs (c++): Completed initialization." << endl;

	unsigned int msg_interval = (100 >= max_iter)?1:100;

	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ 

		if (iter % msg_interval == 0) {
		  cout << "lda_fgs (c++): gibbs iter# " << iter + 1;
		}

		// samples \beta
		for(k = 0; k < num_topics; k++)
			beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));


		for (d = 0; d < num_docs; d++){ // for each document

			vector < unsigned int > word_idx = doc_word_indices[d];

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			theta_samples.col(d) = theta_d;

			// samples z and updates \beta counts
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

			// Note: theta_samples and beta_samples are based on the z from the 
			//       previous iteration   
			
			if (save_z) { Z.col(ss_idx) = z; }
			if (save_beta) { betas.slice(ss_idx) = beta_samples; }
			if (save_theta) { thetas.slice(ss_idx) = theta_samples; }
			if (save_lp) {
			  // Computes log posterior based on (3.4)
			  double logp = calc_log_posterior(theta_samples, beta_samples,
                                         doc_word_indices, doc_lengths, 
                                         word_ids, z, alpha_v, eta);
			  log_posterior(ss_idx) = logp; 
			  if (iter % msg_interval == 0){ cout << " lp: " << logp; }
			}

			ss_idx++; // updates the counter  
			
		} // Handles burn in period


		if (iter % msg_interval == 0) { cout << endl; }

	} // The end of the Gibbs loop
  
	cout << "lda_fgs (c++): Completed Gibbs sampling." << endl;
	cout << "lda_fgs (c++): Number of saved samples - " << ss_idx << endl;
  
	return List::create(
	  Named("thetas") = wrap(thetas),
    Named("betas") = wrap(betas),
		Named("Z") = wrap(Z),
		Named("lp") = wrap(log_posterior)
	);
    
}


/**
 *  The LDA full Gibbs sampler (FGS): 
 *  
 *    It's functionally same as lda_fgs(). It differs only in the input format  
 *    of the corpus, which is based on the Blei corpus format. See: 
 *    http://www.cs.princeton.edu/~blei/lda-c/readme.txt
 *
 *   Arguments:
 * 		num_topics_              - number of topics
 * 		vocab_size_              - vocabulary size
 * 		doc_lengths_             - number of words in each document as a vector
 *    docs_                    - corpus documents (from the Blei corpus)
 * 		topic_assignments_       - initial topic assignment of each word in each corpus document
 * 		alpha_v_                 - hyperparameter vector for the document topic Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - max number of Gibbs iterations to be perfomed
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are stored
 * 		save_z_                  - save the sampled z (N x 1 vector) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_beta_         	     - save the sampled beta (K x V matrix) for each  
 * 		                           iteration, values: {1, 0} 
 * 		save_theta_              - save the sampled theta (K x D matrix) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_lp_                 - compute and save the log posterior of the LDA 
 * 		                           model
 *
 * 	Returns:
 * 		thetas                   - the sampled thetas after the burn in period
 * 		betas                    - the sampled betas after the burn in period
 * 		Z                        - the sampled word topic assignments after the burn in period
 * 		lp                       - the log posterior of the LDA model after ignoring the normalizing constants
 *
 */
RcppExport SEXP lda_fgs_blei_corpus(SEXP num_topics_, SEXP vocab_size_, 
  SEXP doc_lengths_, SEXP docs_, SEXP topic_assignments_, SEXP alpha_v_, 
  SEXP eta_, SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP save_z_, 
  SEXP save_beta_, SEXP save_theta_, SEXP save_lp_) {

	cout << "lda_fgs (c++): init process..." << endl;

	// Variable from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_);
	uvec z = as<uvec>(topic_assignments_);
	vec alpha_v = as<vec>(alpha_v_);

	double eta = as<double>(eta_);
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_lp = as<unsigned int>(save_lp_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = accu(doc_lengths);
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_fgs (c++): the number of saved samples - " << valid_samples << endl;
  cout << "lda_fgs (c++): the number of words in the corpus - " << num_word_instances << endl;

  cube thetas;
  cube betas;
  umat Z; 
  vec log_posterior;
  
  if (save_z) { Z = zeros<umat>(num_word_instances, valid_samples); }
  if (save_beta) { betas = cube(num_topics, vocab_size, valid_samples); }
  if (save_theta) { thetas = cube(num_topics, num_docs, valid_samples); }
  if (save_lp) { log_posterior = zeros<vec>(valid_samples); } 

	mat beta_samples = zeros<mat>(num_topics, vocab_size);
	mat theta_samples = zeros<mat>(num_topics, num_docs);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);

	unsigned int d, i, k, iter, ss_idx = 0, instances = 0, c = 0, word_id, word_count;
	rowvec eta_v = zeros<rowvec>(vocab_size);
	for(k = 0; k < vocab_size; k++)
		eta_v(k) = eta;

  uvec word_ids = zeros<uvec>(num_word_instances);
  vector < vector < unsigned int > > doc_word_indices;
	vector < unsigned int > unique_words, word_idx;


	// Calculates the document word indices

	cout << "Loading documents....";

	for (d = 0; d < num_docs; d++){

    umat document = as<umat>(VECTOR_ELT(docs_, d));
		vector < unsigned int > word_idx;
		
    for (c = 0; c < document.n_cols; c++){
			word_id = document(0,c);
			word_count = document(1,c);
			for (i = 0; i < word_count; i++){
				word_ids(instances) = word_id;
				word_idx.push_back(instances);
				instances++;
			}
		}
		
    doc_word_indices.push_back(word_idx);
	}

	cout << "DONE." << endl;

	// Initilizes beta

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1;
	}

	cout << "lda_fgs (c++): init success..." << endl;

	unsigned int msg_interval = 100;
	if (msg_interval >= max_iter)
		msg_interval = 1;

	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ 

		if (iter % msg_interval == 0)
			cout << "lda_fgs (c++): gibbs iter# " << iter + 1;

		// samples \beta
		for(k = 0; k < num_topics; k++)
			beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));


		for (d = 0; d < num_docs; d++){ // for each document

			vector < unsigned int > word_idx = doc_word_indices[d];

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			theta_samples.col(d) = theta_d;


			// samples z and updates \beta counts
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

		  if (save_z) { Z.col(ss_idx) = z; }
			
			// Note: theta_samples and beta_samples are from old z
			if (save_beta) { betas.slice(ss_idx) = beta_samples; }
			if (save_theta) { thetas.slice(ss_idx) = theta_samples; }

			// Computes log posterior based on (3.4)
			if (save_lp) {
			  double logp = calc_log_posterior(theta_samples, 
                                      beta_samples,
                                      doc_word_indices, 
                                      doc_lengths,
                                      word_ids, 
                                      z, 
                                      alpha_v, 
                                      eta);
			  log_posterior(ss_idx) = logp; 
			  if (iter % msg_interval == 0){ cout << " lp: " << logp; }
			}

			ss_idx++;
		}


		if (iter % msg_interval == 0) cout << endl;

	} // The end of the Gibbs loop
  
	cout << "lda_fgs (c++): the Gibbs sampling is completed." << endl;
	cout << "lda_fgs (c++): the number of saved samples - " << ss_idx << endl;
  
	return List::create(
		Named("thetas") = wrap(thetas),
		Named("betas") = wrap(betas),
		Named("Z") = wrap(Z),
		Named("lp") = wrap(log_posterior)
	);


}

/**
 *  The Augmented Collapsed Gibbs Sampler (ACGS) of LDA:
 * 
 * 		This Gibbs sampler augments the Collapsed Gibbs Sampler (CGS) of LDA 
 *    (Griffiths and Steyvers 2004) that is a Markov chain on Z, with the 
 *    sampling of \beta and \theta variables giving a chain on (Z, Beta, Theta).
 *
 *	References:
 *		1. Finding scientific topics by Griffiths and Steyvers, 2004
 *		2. LDA collapsed Gibbs sampler implementation by David Newman
 *
 *
 * 	Arguments:
 * 		num_topics_              - number of topics
 * 		vocab_size_              - vocabulary size
 * 		word_ids_                - vocabulary ids of each word in each corpus document
 * 		doc_lengths_             - number of words in each document as a vector
 * 		topic_assignments_       - initial topic assignment of each word in each corpus document
 * 		alpha_v_                 - hyperparameter vector for the document topic Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - max number of Gibbs iterations to be perfomed
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are stored
 * 		save_z_                  - save the sampled z (N x 1 vector) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_beta_         	     - save the sampled beta (K x V matrix) for each  
 * 		                           iteration, values: {1, 0} 
 * 		save_theta_              - save the sampled theta (K x D matrix) for each 
 * 		                           saved iteration, values: {1, 0} 
 * 		save_lp_                 - compute and save the log posterior of the LDA 
 * 		                           model
 *
 * 	Returns:
 * 		thetas                   - the sampled theta's after the burn in period
 * 		betas                    - the sampled beta's after the burn in period
 * 		Z                        - the sampled z's (word topic assignments) after the burn in period
 * 		lp                       - the log posterior of the LDA model after ignoring the normalizing constants
 *
 */
RcppExport SEXP lda_acgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
  SEXP doc_lengths_, SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_, 
  SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP save_z_, SEXP save_beta_, 
  SEXP save_theta_, SEXP save_lp_) {

	// variable declarations

	cout << "lda_acgs (c++): init process..." << endl;

	uvec doc_lengths = as<uvec>(doc_lengths_); // the length of each document
	uvec word_ids = as<uvec>(word_ids_); // word indices
	uvec z = as<uvec>(topic_assignments_); // we get this because, we wanna use a given starting point for Gibbs
	vec alpha_v = as<vec>(alpha_v_); // hyperparameters for the document Dirichlets

	double eta = as<double>(eta_); // hyperparameter for the topic Dirichlets
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_lp = as<unsigned int>(save_lp_);
	unsigned int num_docs = doc_lengths.n_elem; // number of documents in the corpus
	unsigned int num_word_instances = word_ids.n_elem; // total number of words in the corpus
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_acgs (c++): the number of saved samples - " << valid_samples 
       << endl;
  cout << "lda_acgs (c++): the number of words in the corpus - " 
       << num_word_instances << endl;

  cube thetas;
  cube betas;
  umat Z; 
  vec log_posterior;
  
  if (save_z) { Z = zeros<umat>(num_word_instances, valid_samples); }
  if (save_beta) { betas = cube(num_topics, vocab_size, valid_samples); }
  if (save_theta) { thetas = cube(num_topics, num_docs, valid_samples); }
  if (save_lp) { log_posterior = zeros<vec>(valid_samples); } 

	mat beta_counts = zeros<mat>(num_topics, vocab_size);
	mat theta_counts = zeros <mat>(num_topics, num_docs);
	vec topic_counts = zeros<vec> (num_topics);
	uvec doc_ids = zeros<uvec>(num_word_instances);
	unsigned int d, i, k, iter, ss_idx = 0, instances = 0, idx;
  unsigned int wid, did, topic, new_topic; 
	double doc_denom;
	vec prob;
	vector < vector < unsigned int > > doc_word_indices;
	mat beta_samples;
	mat theta_samples; 

	// Gets a random permutation of indices
	// this may improve mixing
	// Reference: Dr. Newman's implementation
	// uvec porder = randperm (num_word_instances);


	// Gets the document index for each word instance
	// Calculates the indices for each word in a document as a vector

	for (d = 0; d < num_docs; d++){
		vector < unsigned int > word_idx;
		for (i = 0; i < doc_lengths(d); i++){
			doc_ids(instances) = d;
			word_idx.push_back(instances);
			instances++;
		}
		doc_word_indices.push_back(word_idx);
	}



	// Initializes beta, theta, and topic counts

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (d = 0; d < num_docs; d++)
		theta_counts.col(d) = alpha_v; // initializes with the smoothing parameter

	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1;
		topic_counts (z(i)) += 1;
		theta_counts (z(i), doc_ids(i)) += 1;
	}

	cout << "lda_acgs (c++): init success..." << endl;

	unsigned int msg_interval = 100;
	if (msg_interval >= max_iter) msg_interval = 1;


	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ // for each Gibbs iteration

		if (iter % msg_interval == 0) 
      cout << "lda_acgs (c++): gibbs iter# " << iter + 1;

		
		// Augmenting the collapsed Gibbs sampler chain.
		// It's named as ACGS chain.
		if ((iter >= burn_in) && (iter % spacing == 0)){ // handles the burn in period
		  if (save_beta) { 
		    beta_samples = zeros<mat>(num_topics, vocab_size);
		    for(k = 0; k < num_topics; k++)
		      beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));
		    betas.slice(ss_idx) = beta_samples; 
		  }
		  if (save_theta) { 
		    theta_samples = zeros <mat>(num_topics, num_docs);
		    for (d = 0; d < num_docs; d++) // for each document
		      theta_samples.col(d) = sample_dirichlet(num_topics, theta_counts.col(d));
		    thetas.slice(ss_idx) = theta_samples;
		  }
		}
		
		
		
		for (i = 0; i < num_word_instances; i++){ // for each word instance

			idx = i; //  porder(i); // permutation
			wid = word_ids(idx); // word index
			did = doc_ids(idx); // document index
			topic = z(idx); // old topic
			prob = zeros <vec> (num_topics); // init. probability vector

			// decrements the counts by one, to ignore the current sampling word

			beta_counts(topic, wid)--;
			theta_counts(topic, did)--;
			topic_counts(topic)--;

			doc_denom = doc_lengths(did) - 1 + accu(alpha_v); // a constant for a term

			for (k = 0; k < num_topics; k++){ // for each topic compute P(z_i == j | z_{-i}, w)
				prob(k) = (theta_counts(k, did) / doc_denom) *  (beta_counts(k, wid) /  
                  (topic_counts(k) + vocab_size * eta));
			}
			new_topic = sample_multinomial(prob); // new topic

			// increments the counts by one

			beta_counts(new_topic, wid)++;
			theta_counts(new_topic, did)++;
			topic_counts(new_topic)++;
			z(idx) = new_topic;

		} // end of the word topic sampling loop

		if ((iter >= burn_in) && (iter % spacing == 0)){ // handles the burn in period

		  if (save_z) { Z.col(ss_idx) = z; }

		  // Computes log posterior based on (3.4)
		  if (save_beta && save_theta && save_lp) {
		    double logp = calc_log_posterior(theta_samples, beta_samples, doc_word_indices, 
                                         doc_lengths, word_ids, z, alpha_v, eta);
		    log_posterior(ss_idx) = logp; 
		    if (iter % msg_interval == 0){ cout << " lp: " << logp; }
		  }
		  
		  ss_idx++;

		}

		if (iter % msg_interval == 0) cout << endl;

	} // The end of the Gibbs loop

	cout << "lda_acgs (c++): the Gibbs sampling completed." << endl;
	cout << "lda_acgs (c++): the number of saved samples - " << ss_idx << endl;

	return List::create(
		Named("thetas") = wrap(thetas),
		Named("betas") = wrap(betas),
		Named("Z") = wrap(Z),
		Named("lp") = wrap(log_posterior)
	);

}

