#include "ldamcmc.h"

/**
 * Samples an integer from [0, K) uniformly at random
 * 
 * Arguments: 
 * 		K - the upper interval 
 * Returns:
 * 		the sampled integer  
 */
unsigned int sample_uniform_int(unsigned int K){
	return (unsigned int) (runif(1)(0) * (double)K); // To speedup
}

/**
 * Samples from a given multimomial probability vector
 *
 * Arguments:
 * 		theta - the Multinomial probability vector
 * Returns:
 * 		t - the sampled index
 *
 */
unsigned int sample_multinomial (vec theta) {

	unsigned int t = 0;
	double total_prob = accu(theta);
	double u = runif(1)(0) * total_prob;
	double cumulative_prob = theta(0);

	while(u > cumulative_prob){
		t++;
		cumulative_prob += theta(t);
	}

	return t;

}

/**
 * Samples from a Dirichlet distribution given a set of hyperparameters
 * 
 * Aruguments:
 * 		num_elements - the dimentionality of the Dirichlet distribution 
 * 		alpha - the hyperparameter vector which is in the column vector format  
 * Returns: 
 * 		the Dirichlet sample in the column vector format   
 */
vec sample_dirichlet(unsigned int num_elements, vec alpha){

	vec dirichlet_sample = zeros<vec>(num_elements);

	for ( register unsigned int i = 0; i < num_elements; i++ )
		dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0);

	dirichlet_sample /= accu(dirichlet_sample);

	return dirichlet_sample;

}

/**
 * Samples from a Dirichlet distribution given a set of hyperparameters
 * 
 * Aruguments:
 * 		num_elements - the dimentionality of the Dirichlet distribution 
 * 		alpha - the hyperparameter vector which is in the row vector format  
 * Returns: 
 * 		the Dirichlet sample in the row vector format   
 */
rowvec sample_dirichlet_row_vec (unsigned int num_elements, rowvec alpha){

	rowvec dirichlet_sample = zeros<rowvec>(num_elements);

	for ( register unsigned int i = 0; i < num_elements; i++ )
		dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0);

	dirichlet_sample /= accu(dirichlet_sample);

	return dirichlet_sample;

}

/**
 * Samples random permutations for a given count
 *
 * Arguments:
 * 		n - the number of samples
 * Return:
 * 		order - a vector of indices that represents
 * 				the permutations of numbers in [1, n]
 **/
uvec randperm(unsigned int n) {
	uvec order = zeros<uvec>(n);
	unsigned int k, nn, takeanumber, temp;
	for (k=0; k<n; k++) order(k) = k;
	nn = n;
	for (k=0; k<n; k++) {
		takeanumber = sample_uniform_int(nn); // take a number between 0 and nn-1
		temp = order(nn-1);
		order(nn-1) = order(takeanumber);
		order(takeanumber) = temp;
		nn--;
	}
	return order;
}


vec log_gamma_vec(vec x_vec){

	vec lgamma_vec = zeros<vec>(x_vec.n_elem);

	for (unsigned int i = 0; i < x_vec.n_elem; i++)
		lgamma_vec(i) = lgamma(x_vec(i));

	return lgamma_vec;

}

vec gamma_col_vec(vec x_vec){
  // It took 2hrs of my time in the April 19, 2014 morning to make this function
  // work. The main issue was with accessing the R gamma function from the 
  // RcppArmadillo namespace. See 
  // http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html 
  // gamma(as<NumericVector>(wrap(x_vec))) is another option, but it seems to be
  // slow. See 
  // http://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector

  vec gamma_vec = zeros<vec>(x_vec.n_elem);
  
	for (unsigned int i = 0; i < x_vec.n_elem; i++)
		gamma_vec(i) = Rf_gammafn(x_vec(i));

	return gamma_vec;
}


/***
 * Computes log posterior based on (3.4)
 */
double calc_log_posterior(mat prior_theta, mat prior_beta,
	vector < vector < unsigned int > > doc_word_indices, uvec doc_lengths,
	uvec word_ids, uvec z, vec alpha_v, double eta){

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
 * 		store_dirichlet_         - [1]: saves the sampled (Beta, Theta, Z) from the Markov chain
 * 							                 [0]: saves the sampled Z from the Markov chain
 *
 * 	Returns:
 * 		thetas                   - sampled thetas after the burn in period
 * 		betas                    - sampled betas after the burn in period
 * 		Z                        - sampled word topic assignments after the burn in period
 * 		lp                       - log posterior of the LDA model after ignoring the normalizing constants
 *
 */
RcppExport SEXP lda_fgs(SEXP num_topics_, SEXP vocab_size_, SEXP word_ids_, 
  SEXP doc_lengths_, SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_, 
  SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

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
	unsigned int store_dirichlet = as<unsigned int>(store_dirichlet_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = word_ids.n_elem;
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_fgs (c++): the number of saved samples - " << valid_samples 
       << endl;
  cout << "lda_fgs (c++): the number of words in the corpus - " 
       << num_word_instances << endl;

	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_posterior = zeros<vec>(valid_samples);

	mat prior_beta_samples = zeros<mat>(num_topics, vocab_size);
	mat prior_beta_counts = zeros<mat>(num_topics, vocab_size);
	mat prior_theta_samples = zeros<mat>(num_topics, num_docs);
	mat prior_theta_counts = zeros <mat>(num_topics, num_docs);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);

	vector < vector < unsigned int > > doc_word_indices;
	unsigned int d, i, k, iter, ss_idx = 0, instances = 0;
	rowvec eta_v = zeros<rowvec>(vocab_size);
	for(k = 0; k < vocab_size; k++)
		eta_v(k) = eta;

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

	cout << "lda_fgs (c++): init success..." << endl;

	unsigned int msg_interval = 100;
	if (msg_interval >= max_iter)
		msg_interval = 1;

	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ 

		if (iter % msg_interval == 0) 
      cout << "lda_fgs (c++): gibbs iter# " << iter + 1;

		// samples \beta
		prior_beta_counts = beta_counts; // this is used for log marginal posterior
		for(k = 0; k < num_topics; k++)
			prior_beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, 
                                                           beta_counts.row(k));


		for (d = 0; d < num_docs; d++){ // for each document

			vector < unsigned int > word_idx = doc_word_indices[d];

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			prior_theta_samples.col(d) = theta_d;
			prior_theta_counts.col(d) = partition_counts;


			// samples z and updates \beta counts
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % prior_beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

			// Note: prior_theta_samples and prior_beta_samples are from old z

			Z.col(ss_idx) = z;

			if (store_dirichlet == 1){
				thetas.slice(ss_idx) = prior_theta_samples;
				betas.slice(ss_idx) = prior_beta_samples;
			}



			// Computes log posterior based on (3.4)

			log_posterior(ss_idx) = calc_log_posterior(prior_theta_samples, 
                                                 prior_beta_samples,
                                        				 doc_word_indices, doc_lengths,
                                        				 word_ids, z, alpha_v, eta);


			if (iter % msg_interval == 0) cout << " lp: " << log_posterior(ss_idx);

			ss_idx++;
		}


		if (iter % msg_interval == 0) cout << endl;

	} // The end of the Gibbs loop
  
	cout << "lda_fgs (c++): the Gibbs sampling is completed." << endl;
	cout << "lda_fgs (c++): the number of saved samples - " << ss_idx << endl;
  
	return List::create(Named("thetas") = wrap(thetas),
    Named("betas") = wrap(betas),
		Named("Z") = wrap(Z),
		Named("lp") = wrap(log_posterior));
    
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
 * 		store_dirichlet_         - [1]: saves the sampled (Beta, Theta, Z) from the Markov chain
 * 							                 [0]: saves the sampled Z from the Markov chain
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
  SEXP eta_, SEXP max_iter_, SEXP burn_in_, SEXP spacing_, 
  SEXP store_dirichlet_) {

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
	unsigned int store_dirichlet = as<unsigned int>(store_dirichlet_);

	// Function variables 

	unsigned int num_docs = doc_lengths.n_elem;
	unsigned int num_word_instances = accu(doc_lengths);
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_fgs (c++): the number of saved samples - " << valid_samples << endl;
  cout << "lda_fgs (c++): the number of words in the corpus - " << num_word_instances << endl;

	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_posterior = zeros<vec>(valid_samples);

	mat prior_beta_samples = zeros<mat>(num_topics, vocab_size);
	mat prior_beta_counts = zeros<mat>(num_topics, vocab_size);
	mat prior_theta_samples = zeros<mat>(num_topics, num_docs);
	mat prior_theta_counts = zeros <mat>(num_topics, num_docs);
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
		prior_beta_counts = beta_counts; // this is used for log marginal posterior
		for(k = 0; k < num_topics; k++)
			prior_beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));


		for (d = 0; d < num_docs; d++){ // for each document

			vector < unsigned int > word_idx = doc_word_indices[d];

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			prior_theta_samples.col(d) = theta_d;
			prior_theta_counts.col(d) = partition_counts;


			// samples z and updates \beta counts
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % prior_beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period

			// Note: prior_theta_samples and prior_beta_samples are from old z

			Z.col(ss_idx) = z;

			if (store_dirichlet == 1){
				thetas.slice(ss_idx) = prior_theta_samples;
				betas.slice(ss_idx) = prior_beta_samples;
			}

			// Computes log posterior based on (3.4)

			log_posterior(ss_idx) = calc_log_posterior(prior_theta_samples, 
                                                 prior_beta_samples,
					                                       doc_word_indices, doc_lengths,
					                                       word_ids, z, alpha_v, eta);


			if (iter % msg_interval == 0) cout << " lp: " << log_posterior(ss_idx);

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
		Named("lp") = wrap(log_posterior));


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
 * 		store_dirichlet_         - [1]: saves the ACGS chain
 * 							                 [0]: saves the CGS chain
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
  SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

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
	unsigned int store_dirichlet = as<unsigned int>(store_dirichlet_);
	unsigned int num_docs = doc_lengths.n_elem; // number of documents in the corpus
	unsigned int num_word_instances = word_ids.n_elem; // total number of words in the corpus
	unsigned int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  cout << "lda_acgs (c++): the number of saved samples - " << valid_samples 
       << endl;
  cout << "lda_acgs (c++): the number of words in the corpus - " 
       << num_word_instances << endl;

	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_posterior = zeros<vec>(valid_samples);

	mat beta_counts = zeros<mat>(num_topics, vocab_size);
	mat theta_counts = zeros <mat>(num_topics, num_docs);
	vec topic_counts = zeros<vec> (num_topics);
	uvec doc_ids = zeros<uvec>(num_word_instances);
	unsigned int d, i, k, iter, ss_idx = 0, instances = 0, idx;
  unsigned int wid, did, topic, new_topic; 
	double doc_denom;
	vec prob;
	vector < vector < unsigned int > > doc_word_indices;

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

			Z.col(ss_idx) = z;

			if (store_dirichlet == 1){

				// Augmenting the collapsed Gibbs sampler chain.
				// It's named as ACGS chain.

				mat beta = zeros<mat>(num_topics, vocab_size);
				mat theta = zeros <mat>(num_topics, num_docs);
				for(k = 0; k < num_topics; k++)
					beta.row(k) = sample_dirichlet_row_vec(vocab_size, 
                                                 prior_beta_counts.row(k));
				for (d = 0; d < num_docs; d++) // for each document
					theta.col(d) = sample_dirichlet(num_topics, prior_theta_counts.col(d));
				thetas.slice(ss_idx) = theta;
				betas.slice(ss_idx) = beta;


				// Computes log posterior based on (3.4)

				log_posterior(ss_idx) = calc_log_posterior(theta, beta, 
                                                   doc_word_indices, 
                                                   doc_lengths, word_ids, z,
						                                       alpha_v, eta);

				if (iter % msg_interval == 0) cout << " lp: " << log_posterior(ss_idx);
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
		Named("lp") = wrap(log_posterior));

}

