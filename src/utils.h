///////////////////////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////////////////////

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <math.h>
# include <assert.h>

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

extern
  unsigned int sample_uniform_int (unsigned int K);


extern
  unsigned int sample_multinomial (arma::vec theta);

extern
  arma::vec sample_dirichlet (
      unsigned int num_elements,
      arma::vec alpha
    );

extern
  arma::rowvec sample_dirichlet_row_vec (
      unsigned int num_elements,
      arma::rowvec alpha
    );

extern
  arma::uvec randperm (unsigned int n);

extern
  arma::vec log_gamma_vec (arma::vec x_vec);

extern
  arma::rowvec log_gamma_rowvec(arma::rowvec x_vec);

extern
  arma::vec digamma_vec (arma::vec x_vec);

extern
  arma::rowvec digamma_rowvec(arma::rowvec x_vec);

extern
  arma::vec trigamma_vec (arma::vec x_vec);

extern
  arma::vec tetragamma_vec (arma::vec x_vec);

extern
  arma::vec gamma_col_vec (arma::vec x_vec);

extern
  double sample_antoniak(unsigned int N, double alpha);
