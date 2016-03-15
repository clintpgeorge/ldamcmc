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
  ) {
  
  cout << "lda_fgs (c++): initialization of variables ..." << endl;
  
  // Variable from the R interface  
  
  uvec doc_lengths = as<uvec>(doc_lengths_);
  uvec word_ids = as<uvec>(word_ids_);
  uvec init_z = as<uvec>(topic_assignments_);
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

  uvec z; 
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
    
    // Initilizes beta and z
    
    beta_counts.fill(eta); // initializes with the smoothing parameter
    z = uvec(init_z); // initializes the \z vector 
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
  ) {

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
RcppExport SEXP lda_fgs_blei_corpus(
    SEXP num_topics_, 
    SEXP vocab_size_, 
    SEXP doc_lengths_, 
    SEXP docs_, 
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
  ) {

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
  ) {

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




// =============================================================================
// Serial Tempering 
// =============================================================================

/**
 *  The LDA full Gibbs sampler with Serial Tempering and Tuning  
 *
 *   Arguments:
 *   	num_topics_              - number of topics for the corpus 
 * 		vocab_size_              - vocabulary size of the corpus 
 * 		word_ids_                - a vector of vocabulary ids of word instances in 
 *                               the corpus \eqn{(w_{1,1}, ..., w_{1,n_1}, 
 *                               w_{2,1}, ..., w_{2,n_2}, ..., w_{D,1}, ..., 
 *                               w_{D,n_D})}
 * 		doc_lengths_             - a vector of word counts of documents in the corpus 
 * 		topic_assignments_       - a vector of initial topic assignments for each 
 *                               word in each corpus document
 *    h_grid_                  - a 2-dimensional grid of hyperparameters 
 *                               \eqn{h = (\eta, \alpha)}. It is a 2 x G matrix, 
 *                               where G is the number of grid points and the 
 *                               first row is for \eqn{\alpha} values and the 
 *                               second row is for \eqn{\eta} values  
 *    st_grid_                 - a 2-dimensional grid of hyperparameters 
 *                               \eqn{h = (\eta, \alpha)}. It is a 2 x G matrix, 
 *                               where G is the number of grid points and the 
 *                               first row is for \eqn{\alpha} values and the 
 *                               second row is for \eqn{\eta} values. This a 
 *                               subgrid on h_grid_ that is used for Serial 
 *                               Tempering   
 *    st_grid_nbrs_            - the neighbor indices, from [0, G-1], of each 
 *                               helper grid point
 *   	init_st_grid_index_      - the index of the helper h grid, from [1, G], of 
 *                               the initial hyperparameter \eqn{h = (\eta, 
 *                               \alpha)} 
 *    init_st_grid_zetas_      - initial guess for normalization constants 
 *    max_iter_                - maximum number of Gibbs iterations 
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are saved 
 *    tuning_iter_             - number of tuning iterations 
 * 		save_z_                  - save the sampled Z for each saved iteration, 
 *                               values: {1, 0} 
 * 		save_beta_         	     - save the sampled Beta for each saved iteration, 
 *                               values: {1, 0} 
 * 		save_theta_              - save the sampled Theta for each saved 
 *                               iteration, values: {1, 0} 
 * 		save_st_grid_index_      - save the accepted h's index based on the 
 *                               Metropolis Hastings jump criterion for each 
 *                               saved iteration, values: {1, 0}   
 * 		save_lp_                 - compute and save the log posterior of the LDA 
 *                               model, values: {1, 0} 
 *    save_hat_ratios_         - save h-grid hat ratios, values: {1, 0}  
 *    save_tilde_ratios_       - save h-grid tilde ratios, values: {1, 0}  
 *    verbose_                 - values: {1, 0}   
 *    max_iter_final_          - maximum number of Gibbs iterations for the 
 *                               final run 
 *
 *  Last modified on: 
 *    January 28, 2016 
 */
 
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
  ) {

	cout << endl << "lda_fgs_st (c++): Initializing variables...";

	// Variables from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_);
	uvec word_ids = as<uvec>(word_ids_);
	uvec z = as<uvec>(topic_assignments_);
	mat h_grid = as<mat>(h_grid_);
	mat st_grid = as<mat>(st_grid_);
	vec zetas = as<vec>(init_st_grid_zetas_);
	unsigned int init_st_grid_index = as<unsigned int>(init_st_grid_index_);
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int tuning_iter = as<unsigned int>(tuning_iter_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_st_grid_index = as<unsigned int>(save_st_grid_index_);
	unsigned int save_lp = as<unsigned int>(save_lp_);
	unsigned int save_hat_ratios = as<unsigned int>(save_hat_ratios_);
	unsigned int save_tilde_ratios = as<unsigned int>(save_tilde_ratios_);
	unsigned int verbose = as<unsigned int>(verbose_);
	unsigned int max_iter_final = as<unsigned int>(max_iter_final_);
  
	// Local variables 
	unsigned int d, i, k, g, iter, titer, instances;
	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = word_ids.n_elem;
	// unsigned int num_saved_samples = (save_hat_ratios
  //                                   ? ceil((max_iter-burn_in)*3/(double)spacing)
  //                                  : ceil((max_iter-burn_in)/(double)spacing));
  unsigned int num_saved_samples = ceil((max_iter_final-burn_in)/(double)spacing);
	cube thetas;
	cube betas;
	umat Z; 
  uvec st_grid_index;  
  vec log_posterior;
  vec m_hat;
  mat m_hat_ratios; 
  vec m_tilde;
  mat m_tilde_ratios; 
  vec st_grid_m_hat; 

	if (save_z) { Z = zeros<umat>(num_word_instances, num_saved_samples); }
	if (save_beta) { betas = cube(num_topics, vocab_size, num_saved_samples); }
	if (save_theta) { thetas = cube(num_topics, num_docs, num_saved_samples); }
  if (save_st_grid_index) { st_grid_index = zeros<uvec>(num_saved_samples); }
  if (save_lp) { log_posterior = zeros<vec>(num_saved_samples); } 
  if (save_hat_ratios){ m_hat_ratios = zeros<mat>(h_grid.n_cols, num_saved_samples); }
  if (save_tilde_ratios){ m_tilde_ratios = zeros<mat>(h_grid.n_cols, num_saved_samples);}

	mat beta_samples = zeros<mat>(num_topics, vocab_size);
	mat theta_samples = zeros<mat>(num_topics, num_docs);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);
  


  // Gets the neighbours of each point in st.grid, from the R list object 
  vector < vector < unsigned int > > st_grid_nbr_indices;
  for (g = 0; g < st_grid.n_cols; g++){
    vector < unsigned int > nbr_indices; 
    rowvec neighbors = as<rowvec>(VECTOR_ELT(st_grid_nbrs_, g));
    for (i = 0; i < neighbors.n_elem; i++){
      nbr_indices.push_back(neighbors(i));
    }
    st_grid_nbr_indices.push_back(nbr_indices);
  }


	// Calculates the indices for words in every document, as a vector
  vector < vector < unsigned int > > doc_word_indices;
  instances = 0;
	for (d = 0; d < num_docs; d++){
		vector < unsigned int > word_idx;
		for (i = 0; i < doc_lengths(d); i++){
			word_idx.push_back(instances);
			instances++;
		}
		doc_word_indices.push_back(word_idx);
	}

	// Initializes beta counts 
	for (i = 0; i < num_word_instances; i++){ 
		beta_counts(z(i), word_ids(i)) += 1.;
	}

	cout << " SUCCESS." << endl << endl;



  // ***************************************************************************
  // BEGIN: Serial Tempering 
  // ***************************************************************************
  
  vector <unsigned int> curr_h_nbr_indices; 
  vector <unsigned int> prop_h_nbr_indices; 
  unsigned int curr_h_index = init_st_grid_index; 
  unsigned int prop_h_index; 
  unsigned int ss_idx;
  unsigned int num_gibbs_iter;
  double pd_ratio, prior_ratio, zeta_ratio, r, alpha, eta, num_accept; 
  // double MIN_PRIOR_RATIO = 1e+0;
  double M_SMOOTHING_CONST = 1e-20;
  double st_grid_nc; // stores st-grid normalizing constant for an iteration  
  vec h_grid_mh; // stores h-grid \hat{M}(h) ratios for an iteration 
  vec h_grid_mt; // stores h-grid \tilde{M}(h) ratios for an iteration
  vec st_grid_mh; // stores st-grid \hat{M}(h) ratios for an iteration  
  vec st_grid_occupancy;
  mat st_grid_occupancies = zeros<mat>(st_grid.n_cols, tuning_iter); 
  mat st_grid_zetas = zeros<mat>(zetas.n_elem, tuning_iter);
  vec st_grid_pr = zeros<vec>(st_grid.n_cols); // stores prior ratios for st-grid 
  vec h_grid_pr = zeros<vec>(h_grid.n_cols); // stores prior ratios for h-grid 
  int mh_inf_count;
    
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

    // We triple the number of Gibbs iterations for the last tuning iteration 
    // Added on January 23, 2015 
    // num_gibbs_iter = ((((titer+1) == tuning_iter) && save_hat_ratios)
    //                   ? ((max_iter-burn_in)*3 + burn_in) : max_iter);
    // num_gibbs_iter = max_iter; 
    
    // Added on January 28, 2016 
    num_gibbs_iter = (((titer+1) == tuning_iter) ? max_iter_final : max_iter);
    
    cout << "lda_fgs_st (c++): Gibbs sampling (iterations = " << num_gibbs_iter << ")" << endl;
    
    // *************************************************************************
    // BEGIN: Gibbs Sampling Loop 
    // *************************************************************************
  	for (iter = 0; iter < num_gibbs_iter; iter++){ 
  
  		if (verbose) cout << "lda_fgs_st (c++): iter# " << iter + 1;
      
      alpha = st_grid(0, curr_h_index); // alpha from MH 
      eta = st_grid(1, curr_h_index); // alpha from MH 
      
  		// samples \beta
  		for(k = 0; k < num_topics; k++)
  			beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k) + eta);
  
  		for (d = 0; d < num_docs; d++){ // for each document
  
  			vector < unsigned int > word_idx = doc_word_indices[d];
  
  			// samples \theta
        // TODO: need to check whether removing partition_counts speeds up  
  			vec partition_counts = zeros<vec>(num_topics); 
  			for (i = 0; i < doc_lengths(d); i++) 
  				partition_counts(z(word_idx[i])) += 1.;
  			vec theta_d = sample_dirichlet(num_topics, partition_counts + alpha);
  			theta_samples.col(d) = theta_d;

  
  			// samples z and updates \beta counts
        // TODO: need to check why there are 3 for loops 
  			for(i = 0; i < doc_lengths(d); i++)
  				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1.; // excludes document d's word-topic counts
  			for (i = 0; i < doc_lengths(d); i++)
  				z(word_idx[i]) = sample_multinomial(theta_d % beta_samples.col(word_ids(word_idx[i])));
  			for(i = 0; i < doc_lengths(d); i++)
  				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1.; // includes document d's word-topic counts
  
  		}
  
  		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

        // *********************************************************************
        // DO NOT CHANGE THE ORDER OF THIS CODE BLOCK 
        // Dependency: h_grid_mh calculation
        // *********************************************************************

        // Calculates \nu_h(\psi) / \nu_{h_(i)}(\psi) for st-grid  
        // Note: We compute ratios of priors instead of actual prior because  
        // sometimes prior may become too large to handle. 
        st_grid_nc = 0.0; // normalizing constant 
        for (g = 0; g < st_grid.n_cols; g++){
          prior_ratio = calc_prior_ratio(beta_samples, 
                                         theta_samples, 
                                         st_grid(0, g), // alpha 
                                         st_grid(1, g), // eta 
                                         alpha, // st alpha 
                                         eta); // st eta  
          // prior_ratio = is_finite(prior_ratio)?prior_ratio:MIN_PRIOR_RATIO; // This is a HACK!!
          st_grid_pr(g) = prior_ratio;
          st_grid_nc += (prior_ratio / ((double)st_grid.n_cols * zetas(g)));
        }

        if (verbose) cout << " st_grid_nc: " << st_grid_nc;

        // Computes the st-grid \hat{M}(h) as an online average 
        assert(st_grid_nc > 0);
        st_grid_mh = st_grid_pr / st_grid_nc;
        st_grid_m_hat = ((ss_idx * st_grid_m_hat + st_grid_mh) / (ss_idx + 1.)); 
        
        

        // *********************************************************************
        // DO NOT CHANGE THE ORDER OF THIS CODE BLOCK 
        // *********************************************************************
  
        if ((titer + 1) == tuning_iter){ // THE LAST TUNING ITERARION 

          if (save_z) { Z.col(ss_idx) = z; }
          if (save_beta) { betas.slice(ss_idx) = beta_samples; }
          if (save_theta) { thetas.slice(ss_idx) = theta_samples; }
          if (save_st_grid_index) { st_grid_index(ss_idx) = curr_h_index; }
          

          // Calculates \nu_h(\psi) / \nu_{h_(i)}(\psi) for h-grid  
          // Computes the prior ratios for every h in the grid given the h = 
          // (eta, alpha) in ST 
          int inf_count = 0;
          for (g = 0; g < h_grid.n_cols; g++){
            prior_ratio = calc_prior_ratio(beta_samples, 
                                           theta_samples,
                                           h_grid(0, g), 
                                           h_grid(1, g), 
                                           alpha, 
                                           eta);
            // prior_ratio = is_finite(prior_ratio)?prior_ratio:MIN_PRIOR_RATIO; // This is a HACK!!
            if (!is_finite(prior_ratio)) inf_count++;
            h_grid_pr(g) = prior_ratio;
          }
          
          if (verbose) { cout << " count(NaNs): " << inf_count; }
          
          // Computes the h-grid \hat{M}(h) as an online average 
          h_grid_mh =  h_grid_pr / st_grid_nc; //st_grid_nc is *from* st-grid 
          m_hat = ((ss_idx * m_hat + h_grid_mh) / (ss_idx + 1.));  

          // Computes the h-grid \tilde{M}(h) as an online average 
          h_grid_mt = h_grid_pr * zetas(curr_h_index);     
          m_tilde = ((ss_idx * m_tilde + h_grid_mt) / (ss_idx + 1.));  
		  
          // Saves h-grid \hat{M} and \tilde{M} ratios:
          if (save_hat_ratios){ 
            m_hat_ratios.col(ss_idx) = h_grid_mh; 
          }
          if (save_tilde_ratios){
            m_tilde_ratios.col(ss_idx) = h_grid_mt;
          } 
		  
          // Computes log posterior based on (3.4)
          if (save_lp) {
            vec alpha_v = zeros<vec>(num_topics); 
            alpha_v.fill(alpha);
            double logp = calc_log_posterior(theta_samples, beta_samples,
                                             doc_word_indices, doc_lengths, 
                                             word_ids, z, alpha_v, eta);
            log_posterior(ss_idx) = logp; 
            if (verbose){ cout << " lp: " << logp; }
          }
        
        }  // THE LAST TUNING ITERARION 


        ss_idx++; // updates the counter  
  		
      } // Handles burn in period
  


      // ***********************************************************************
      // BEGIN: Metropolis Hastings/Serial Tempering Jump
      // ***********************************************************************
       
      st_grid_occupancy(curr_h_index) += 1.; // updates helper-grid occupancy
      
      curr_h_nbr_indices = st_grid_nbr_indices[curr_h_index];
      prop_h_index = curr_h_nbr_indices[sample_uniform_int(curr_h_nbr_indices.size())];
      prop_h_nbr_indices = st_grid_nbr_indices[prop_h_index];
      
      // Computing the acceptance ratio 
      pd_ratio = ((double)curr_h_nbr_indices.size()
                  / (double)prop_h_nbr_indices.size()); 
      prior_ratio = calc_prior_ratio(beta_samples, 
                                     theta_samples, 
                                     st_grid(0, prop_h_index), // proposed alpha 
                                     st_grid(1, prop_h_index), // proposed eta
                                     alpha, // current alpha 
                                     eta); // current eta 
      // TODO: sometimes the prior ratio can be 0 or infinity
      // How do we handle it? Not sure yet!!!! 
      if (is_finite(prior_ratio)){
        zeta_ratio = zetas(curr_h_index) / zetas(prop_h_index); 
        r = min(1.0000, pd_ratio * prior_ratio * zeta_ratio);
        if (verbose) {
          printf(" (r=%.3f,", r);
          printf(" pdr=%.3f,", pd_ratio); // proposal density ratio 
          printf(" ppr=%.3f,", zeta_ratio); // the ratio of normalizing constants  
          printf(" pr=%.3f)", prior_ratio); // Bayes factor ratio 
        }
        if (runif(1)(0) < r){ // acceptance 
          curr_h_index = prop_h_index; // updates the current h index  
          cout << (verbose ? " A" : "+");
          num_accept += 1.; 
        }
        else {
          cout << (verbose ? " R" : "-");
        }
      }
      else { 
        // Conditional acceptance --- This is a HACK
        curr_h_index = prop_h_index; // updates the current h index  
        cout << (verbose ? " A*" : "*");
        num_accept += 1.;       
        mh_inf_count++; 
      }
      
      // ***********************************************************************
      // END: Metropolis Hastings/Serial Tempering Jump
      // ***********************************************************************
      if (verbose) cout << endl;
  
  	} 
    // *************************************************************************
    // END: Gibbs Sampling Loop 
    // *************************************************************************
    cout << endl;
  	cout << "lda_fgs_st (c++): End of Gibbs sampling" << endl;
  	cout << "lda_fgs_st (c++): Number of saved samples = " << ss_idx << endl;
    cout << "lda_fgs_st (c++): Serial Tempering acceptance (%) = " 
         << 100. * (num_accept / num_gibbs_iter) << endl;
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
    
//  	if (min(st_grid_m_hat) <= 0.){ 
//      // a safe exit from the serial tempering chain 
//      // TODO: need to check whether a smoothing on Bayes factors helps 
//      // regularization? 
//      cout << "lda_fgs_st (c++): iter# " << iter + 1 << endl;
//      cout << "lda_fgs_st (c++): Exiting the serial tempering chain to avoid DBZ." << endl;
//      break; 
//  	}
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
    Named("thetas") = wrap(thetas),
    Named("betas") = wrap(betas),
    Named("Z") = wrap(Z),
    Named("lp") = wrap(log_posterior),    
    Named("m_hat_ratios") = wrap(m_hat_ratios),
    Named("m_hat") = wrap(m_hat),
    Named("m_tilde") = wrap(m_tilde),
    Named("st_grid_index") = wrap(st_grid_index),
    Named("st_grid_zetas") = wrap(st_grid_zetas), 
    Named("st_grid_occupancies") = wrap(st_grid_occupancies),
    Named("st_grid_m_hat") = wrap(st_grid_m_hat),
    Named("m_tilde_ratios") = wrap(m_tilde_ratios)
  );

}

/**
 *  The LDA full gibbs sampler with hyperparameter selection:
 *
 * 	Arguments:
 * 		num_topics_              - number of topics
 * 		vocab_size_              - vocabulary size
 * 		word_ids_                - vocabulary ids of each word in each corpus document
 * 		doc_lengths_             - number of words in each document as a vector
 * 		topic_assignments_       - initial topic assignment of each word in each corpus document
 * 		alpha_                   - hyperparameter for the document topic Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - maximum number of Gibbs iterations to be perfomed
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are saved 
 * 		save_z_                  - save the sampled Z for each saved iteration, values: {1, 0} 
 * 		save_beta_         	     - save the sampled Beta for each saved iteration, values: {1, 0} 
 * 		save_theta_              - save the sampled Theta for each saved iteration, values: {1, 0} 
 * 		save_Bh_                 - save the computed B(h, h_*) for each saved iteration, values: {1, 0}   
 * 		save_lp_                 - compute and save the log posterior of the LDA model
 *
 * 	Returns:
 * 		thetas                   - the sampled thetas after the burn in period
 * 		betas                    - the sampled betas after the burn in period
 * 		Z                        - the sampled word topic assignments after the burn in period
 * 		lp                       - the log posterior of the LDA model after ignoring the normalizing constants
 *
 */
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
  ) {

	cout << "lda_fgs_hs (c++): init process..." << endl;

	// Variables from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_);
	uvec word_ids = as<uvec>(word_ids_);
	uvec z = as<uvec>(topic_assignments_);
	mat h_grid = as<mat>(h_grid_);
	double alpha = as<double>(alpha_);
	double eta = as<double>(eta_);
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_Bh = as<unsigned int>(save_Bh_);
	unsigned int save_lp = as<unsigned int>(save_lp_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = word_ids.n_elem;
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
	cout << "lda_fgs_hs (c++): number of saved samples - " << valid_samples 
			<< endl;
	cout << "lda_fgs_hs (c++): number of words in the corpus - " 
			<< num_word_instances << endl;

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

	vec alpha_v = zeros<vec>(num_topics); 
	alpha_v.fill(alpha);
  
	vector < vector < unsigned int > > doc_word_indices;
	unsigned int d, i, k, g, iter, ss_idx = 0, instances = 0;


	// Calculates the indices for each word in a document as a vector
	for (d = 0; d < num_docs; d++){
		vector < unsigned int > word_idx;
		for (i = 0; i < doc_lengths(d); i++){
			word_idx.push_back(instances);
			instances++;
		}
		doc_word_indices.push_back(word_idx);
	}
 
	uvec h_grid_sort_idx; 
	vec h_grid_mean_bf = zeros<vec>(h_grid.n_cols);  
  vec h_grid_bf = zeros<vec>(h_grid.n_cols);
	mat h_grid_bf_s; 
	if (save_Bh) { h_grid_bf_s = zeros<mat>(h_grid.n_cols, valid_samples); }

	

	// Initilizes beta

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1.;
	}

	cout << "lda_fgs_hs (c++): init success..." << endl;

	unsigned int msg_interval = 100;
	if (msg_interval >= max_iter) { msg_interval = 1; }

	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ 

		if (iter % msg_interval == 0)
			cout << "lda_fgs_hs (c++): gibbs iter# " << iter + 1;

		// samples \beta
		for(k = 0; k < num_topics; k++)
			beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));


		for (d = 0; d < num_docs; d++){ // for each document

			vector < unsigned int > word_idx = doc_word_indices[d];

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1.;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			theta_samples.col(d) = theta_d;


			// samples z and updates \beta counts
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1.; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1.; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

			// Note: theta_samples and beta_samples are from old z
			if (save_z) { Z.col(ss_idx) = z; }
			if (save_beta) { betas.slice(ss_idx) = beta_samples; }
			if (save_theta) { thetas.slice(ss_idx) = theta_samples; }
      
      
			// Computing average of Bayes Factor ratios via incremental average. This
			// is a way to save the space required to store \beta and \theta samples. 
			// See
			// http://math.stackexchange.com/questions/106700/incremental-averageing
			// incremental average. 

      for (g = 0; g < h_grid.n_cols; g++){
        h_grid_bf(g) = calc_prior_ratio(beta_samples, 
                                        theta_samples, 
                                        h_grid(0, g), 
                                        h_grid(1, g), 
                                        alpha, 
                                        eta);
      }
			if (save_Bh){ h_grid_bf_s.col(ss_idx) = h_grid_bf; }			
			h_grid_mean_bf = (ss_idx * h_grid_mean_bf + h_grid_bf) / (ss_idx + 1);  
			h_grid_sort_idx = sort_index(h_grid_mean_bf, "descend");

			// Computes log posterior based on (3.4)

			double logp = calc_log_posterior(theta_samples, beta_samples,
												               doc_word_indices, doc_lengths, word_ids, 
                                       z, alpha_v, eta);
			if (save_lp) { log_posterior(ss_idx) = logp; }

			if (iter % msg_interval == 0){
		        cout << " lp: " << logp; 
		        cout << " h: (" << h_grid(0, h_grid_sort_idx(0)); // shows the max value 
		        cout << ", " << h_grid(1, h_grid_sort_idx(0));
		        cout << ") avg. B(h, h*):" << h_grid_mean_bf(h_grid_sort_idx(0));
			}

			ss_idx++;
		}


		if (iter % msg_interval == 0) cout << endl;

	} // The end of the Gibbs loop
  
	cout << "lda_fgs_hs (c++): End of Gibbs sampling." << endl;
	cout << "lda_fgs_hs (c++): number of saved samples - " << ss_idx << endl;
  
  return List::create(
    Named("thetas") = wrap(thetas),
    Named("betas") = wrap(betas),
    Named("Z") = wrap(Z),
    Named("lp") = wrap(log_posterior), 
    Named("h_grid_mean_bf") = wrap(h_grid_mean_bf), 
    Named("h_grid_sort_idx") = wrap(h_grid_sort_idx), // index starts at 0 
    Named("h_grid_bf_s") = wrap(h_grid_bf_s)
    );

}

/**
 *  The Augmented Collapsed Gibbs Sampler (ACGS) of LDA  with hyperparameter 
 *  selection:
 *    This augments the Collapsed Gibbs Sampler (CGS) of
 * 		LDA (Griffiths and Steyvers 2004), which is a Markov
 * 		chain on Z, with the sampling of \beta and \theta
 * 		variables giving a chain on (Z, Beta, Theta).
 *
 *	References:
 *		1. Finding scientific topics by Griffiths and Steyvers, 2004
 *		2. LDA collapsed Gibbs sampler implementation by David Newman
 *
 *
 * 	Arguments:
 * 		num_topics_              - the number of topics
 * 		vocab_size_              - the vocabulary size
 * 		word_ids_                - the vocabulary ids of each word in each 
 *                               corpus document
 * 		doc_lengths_             - the number of words in each document as a 
 *                               vector
 * 		topic_assignments_       - the initial topic assignment of each word 
 *                               in each corpus document
 * 		alpha_                   - hyperparameter for the document topic 
 *                               Dirichlet
 * 		eta_                     - hyperparameter for the topic Dirichlet
 * 		max_iter_                - maximum number of Gibbs iterations 
 * 		burn_in_                 - burn in period
 * 		spacing_                 - spacing between samples that are saved 
 * 		save_z_                  - save the sample Z for each saved iteration, values: {1, 0} 
 * 		save_beta_         	     - save the sample Beta for each saved iteration, values: {1, 0} 
 * 		save_theta_              - save the sample Theta for each saved iteration, values: {1, 0} 
 * 		save_Bh_                 - save the computed B(h, h_*) for each saved iteration, values: {1, 0}   
 * 		save_lp_                 - compute and save the log posterior of the LDA model
 *
 *
 * 	Returns:
 * 		thetas                   - the sampled theta's after the burn in period
 * 		betas                    - the sampled beta's after the burn in period
 * 		Z                        - the sampled z's (word topic assignments) after the burn in period
 * 		lp                       - the log posterior of the LDA model after ignoring the normalizing constants
 *
 */
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
  ) {

	cout << "lda_acgs_hs (c++): init process..." << endl;

	// Variables from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_); // the length of each document
	uvec word_ids = as<uvec>(word_ids_); // word indices
	uvec z = as<uvec>(topic_assignments_); // the starting point for Gibbs
	mat h_grid = as<mat>(h_grid_);
	double alpha = as<double>(alpha_);
	double eta = as<double>(eta_); // hyperparameter for the topic Dirichlets
	unsigned int num_topics = as<unsigned int>(num_topics_);
	unsigned int vocab_size = as<unsigned int>(vocab_size_);
	unsigned int max_iter = as<unsigned int>(max_iter_);
	unsigned int burn_in = as<unsigned int>(burn_in_);
	unsigned int spacing = as<unsigned int>(spacing_);
	unsigned int save_z = as<unsigned int>(save_z_);
	unsigned int save_beta = as<unsigned int>(save_beta_);
	unsigned int save_theta = as<unsigned int>(save_theta_);
	unsigned int save_Bh = as<unsigned int>(save_Bh_);
	unsigned int save_lp = as<unsigned int>(save_lp_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem; // number of documents in the corpus
	unsigned int num_word_instances = word_ids.n_elem; // total number of words in the corpus
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);

	cout << "lda_acgs_hs (c++): the number of words in the corpus - " 
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
	vec alpha_v = zeros<vec>(num_topics); 
	alpha_v.fill(alpha);

	vec topic_counts = zeros<vec> (num_topics);
	uvec doc_ids = zeros<uvec>(num_word_instances);
	unsigned int d, i, k, g, iter, idx, ss_idx = 0, instances = 0;
	unsigned int wid, did, topic, new_topic;
	double doc_denom;
	vec prob;
	vector < vector < unsigned int > > doc_word_indices;

	// Variables for hyperparameter selection

	uvec h_grid_sort_idx; 
	vec h_grid_mean_bf = zeros<vec>(h_grid.n_cols);
	vec h_grid_bf = zeros<vec>(h_grid.n_cols);
	mat h_grid_bf_s; 	
	if (save_Bh) { h_grid_bf_s = zeros<mat>(h_grid.n_cols, valid_samples); }

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

	beta_counts.fill(eta); // initializes with the smoothing parameter eta
	for (d = 0; d < num_docs; d++)
		theta_counts.col(d) = alpha_v; // initializes with alpha 

	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1.;
		topic_counts (z(i)) += 1.;
		theta_counts (z(i), doc_ids(i)) += 1.;
	}

	cout << "lda_acgs_hs (c++): init success..." << endl;

	unsigned int msg_interval = 100;
	if (msg_interval >= max_iter) { msg_interval = 1; }


	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ // for each Gibbs iteration

		if (iter % msg_interval == 0) 
			cout << "lda_acgs_hs (c++): gibbs iter# " << iter + 1;

		mat prior_beta_counts = beta_counts;
		mat prior_theta_counts = theta_counts;

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

			// a constant for a term
			doc_denom = doc_lengths(did) - 1 + accu(alpha_v); 

			// for each topic compute P(z_i == j | z_{-i}, w)
			for (k = 0; k < num_topics; k++){ 
				prob(k) = ((theta_counts(k, did) / doc_denom) 
							* (beta_counts(k, wid) 
								/ (topic_counts(k) + vocab_size * eta)));
			}
			new_topic = sample_multinomial(prob); // new topic

			// increments the counts by one

			beta_counts(new_topic, wid)++;
			theta_counts(new_topic, did)++;
			topic_counts(new_topic)++;
			z(idx) = new_topic;

		} // end of the word topic sampling loop

		if ((iter >= burn_in) && (iter % spacing == 0)){ // handles the burn in period


			// Augmenting the collapsed Gibbs sampler chain.
			// It's named as ACGS chain.

			mat prior_beta_samples = zeros<mat>(num_topics, vocab_size);
			mat prior_theta_samples = zeros <mat>(num_topics, num_docs);
			for(k = 0; k < num_topics; k++) // for each topic 
				prior_beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, prior_beta_counts.row(k));
			for (d = 0; d < num_docs; d++) // for each document
				prior_theta_samples.col(d) = sample_dirichlet(num_topics, prior_theta_counts.col(d));

			// Saves beta, theta, and z 
			// Note: prior_theta_samples and prior_beta_samples are from old z

			if (save_z) { Z.col(ss_idx) = z; }
			if (save_beta) { betas.slice(ss_idx) = prior_beta_samples; }
			if (save_theta) { thetas.slice(ss_idx) = prior_theta_samples; }
      
      
			// Computing average of Bayes Factor ratios via incremental 
			// average. This is a way to save the space required to store 
			// \beta and \theta samples. See
			// http://math.stackexchange.com/questions/106700/incremental-averageing
			// incremental average. 
			
			for (g = 0; g < h_grid.n_cols; g++){
			  h_grid_bf(g) = calc_prior_ratio(prior_beta_samples, 
               prior_theta_samples, 
               h_grid(0, g), 
               h_grid(1, g), 
               alpha, 
               eta);
			}
			if (save_Bh){ h_grid_bf_s.col(ss_idx) = h_grid_bf; }			
			h_grid_mean_bf = (ss_idx * h_grid_mean_bf + h_grid_bf) / (ss_idx + 1);  
			h_grid_sort_idx = sort_index(h_grid_mean_bf, "descend");
			
			
			// Computes log posterior based on (3.4)

			double logp = calc_log_posterior(prior_theta_samples, prior_beta_samples,
                                        doc_word_indices, doc_lengths, word_ids, 
                                        z, alpha_v, eta);
			if (save_lp) { log_posterior(ss_idx) = logp; }

			if (iter % msg_interval == 0){
        cout << " lp: " << logp; 
        cout << " h: (" << h_grid(0, h_grid_sort_idx(0)); 
        cout << ", " << h_grid(1, h_grid_sort_idx(0));
        cout << ") avg. B(h, h*):"; 
        cout << h_grid_mean_bf(h_grid_sort_idx(0));
			}

			ss_idx++;

		}

		if (iter % msg_interval == 0) cout << endl;

	} // The end of the Gibbs loop

	cout << "lda_acgs_hs (c++): the Gibbs sampling is completed." << endl;
	cout << "lda_acgs_hs (c++): the number of saved samples - " << ss_idx 
	     << endl;

	return List::create(
		Named("thetas") = wrap(thetas),
		Named("betas") = wrap(betas),
		Named("Z") = wrap(Z), // index starts at 0 
		Named("lp") = wrap(log_posterior), 
		Named("h_grid_mean_bf") = wrap(h_grid_mean_bf), 
		Named("h_grid_sort_idx") = wrap(h_grid_sort_idx), // index starts at 0 
		Named("h_grid_bf_s") = wrap(h_grid_bf_s)
		);

}
