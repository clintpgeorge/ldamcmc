///////////////////////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////////////////////

# include "utils.h"




//' Samples from the Antoniak distribution
//'
//' It's done by sampling \eqn{N} Bernoulli variables
//'
//' References:
//'
//'   http://www.jmlr.org/papers/volume10/newman09a/newman09a.pdf
//'
//' @param N Number of samples
//' @param alpha strength parameter
//'
//' @export
//'
//' @family utils
//'
//' @note
//'
//' Created on: May 19, 2016
//'
//' Created by: Clint P. George
//'
// [[Rcpp::export]]
double sample_antoniak(unsigned int N, double alpha){
  vec bs = zeros<vec>(N);
  for (unsigned int l = 0; l < N; l++){
    bs(l) = rbinom(1, 1, (alpha / (alpha + l)))(0);
  }
  return sum(bs);
}

/**
* Samples an integer from [0, K) uniformly at random
*
* Arguments:
* 		K - the upper interval
* Returns:
* 		the sampled integer
*/
unsigned int sample_uniform_int (unsigned int K){
  return (unsigned int) (runif(1)(0) * (double)K); // To speedup
}

//' A speedy sampling from a multimomial distribution
//'
//' @param theta a multinomial probability vector (K x 1 vector)
//'
//' @return returns a class index from [0, K)
//'
//' @note
//' Author: Clint P. George
//'
//' Created on: February 11, 2016
//'
//' @family utils
//'
//' @export
// [[Rcpp::export]]
unsigned int sample_multinomial (arma::vec theta) {

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


//' Samples from a Dirichlet distribution given a hyperparameter
//'
//' @param num_elements the dimention of the Dirichlet distribution
//' @param alpha the hyperparameter vector (a column vector)
//'
//' @return returns a Dirichlet sample (a column vector)
//'
//' @note
//' Author: Clint P. George
//'
//' Created on: 2014
//'
//' @family utils
//'
//' @export
// [[Rcpp::export]]
arma::vec sample_dirichlet(unsigned int num_elements, arma::vec alpha){

  arma::vec dirichlet_sample = arma::zeros<arma::vec>(num_elements);

  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0); // R::rgamma(1, alpha(i));

  dirichlet_sample /= accu(dirichlet_sample);

  return dirichlet_sample;

}

/**
* Samples from a Dirichlet distribution given a hyperparameter
*
* Aruguments:
* 		num_elements - the dimention of the Dirichlet distribution
* 		alpha - the hyperparameter vector (a column vector)
* Returns:
* 		the Dirichlet sample (a column vector)
*/
arma::rowvec sample_dirichlet_row_vec (unsigned int num_elements, arma::rowvec alpha){

  arma::rowvec dirichlet_sample = arma::zeros<arma::rowvec>(num_elements);

  for ( register unsigned int i = 0; i < num_elements; i++ )
    dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0); // R::rgamma(1, alpha(i));

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
arma::uvec randperm(unsigned int n) {
  arma::uvec order = arma::zeros<arma::uvec>(n);
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


arma::vec log_gamma_vec(arma::vec x_vec) {
  arma::vec lgamma_vec = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    lgamma_vec(i) = lgamma(x_vec(i));
  return lgamma_vec;
}

arma::rowvec log_gamma_rowvec(arma::rowvec x_vec) {
  arma::rowvec lgamma_rowvec = arma::zeros<arma::rowvec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    lgamma_rowvec(i) = lgamma(x_vec(i));
  return lgamma_rowvec;
}



arma::vec digamma_vec(arma::vec x_vec) {
  // digamma(wrap()) will do, with comparable performance
  arma::vec ret = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_digamma(x_vec(i));
  return ret;
}

arma::rowvec digamma_rowvec(arma::rowvec x_vec) {
  arma::rowvec ret = arma::zeros<arma::rowvec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_digamma(x_vec(i));
  return ret;
}

arma::vec trigamma_vec(arma::vec x_vec) {
  arma::vec ret = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_trigamma(x_vec(i));
  return ret;
}

arma::vec tetragamma_vec(arma::vec x_vec) {
  arma::vec ret = arma::zeros<arma::vec>(x_vec.n_elem);
  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    ret(i) = Rf_tetragamma(x_vec(i));
  return ret;
}

arma::vec gamma_col_vec(arma::vec x_vec){
  // It took 2hrs of my time in the April 19, 2014 morning to make this function
  // work. The main issue was with accessing the R gamma function from the
  // RcppArmadillo namespace. See
  // http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html
  // gamma(as<NumericVector>(wrap(x_vec))) is another option, but it seems to be
  // slow. See
  // http://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector

  arma::vec gamma_vec = arma::zeros<arma::vec>(x_vec.n_elem);

  for (unsigned int i = 0; i < x_vec.n_elem; i++)
    gamma_vec(i) = Rf_gammafn(x_vec(i));

  return gamma_vec;
}
