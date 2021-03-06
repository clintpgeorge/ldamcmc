% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lda_fgs_ppc}
\alias{lda_fgs_ppc}
\title{LDA: Full Gibbs Sampler with Posterior Predictive Value}
\usage{
lda_fgs_ppc(num_topics, vocab_size, docs_tf, alpha_h, eta_h, max_iter, burn_in,
  spacing, verbose)
}
\arguments{
\item{num_topics}{Number of topics in the corpus}

\item{vocab_size}{Vocabulary size}

\item{docs_tf}{A list of corpus documents read from the Blei corpus using
\code{\link{read_docs}} (term indices starts with 0)}

\item{alpha_h}{Hyperparameter for \eqn{\theta} sampling}

\item{eta_h}{Smoothing parameter for the \eqn{\beta} matrix}

\item{max_iter}{Maximum number of Gibbs iterations to be performed}

\item{burn_in}{Burn-in-period for the Gibbs sampler}

\item{spacing}{Spacing between the stored samples (to reduce correlation)}

\item{verbose}{from {0, 1, 2}}
}
\value{
The Markov chain output as a list of
  \item{lppv}{log posterior predictive values of each document}
  \item{lppc}{averge of log posterior predictive values}
}
\description{
Implements the Full Gibbs sampler for the LDA model---a Markov chain on
\eqn{(\beta, \theta, z)}. The log posterior predictive value is based on 
Zhe Chen (2015)
}
\seealso{
Other MCMC: \code{\link{lda_acgs_st}},
  \code{\link{lda_cgs_em_perplexity}},
  \code{\link{lda_cgs_em}},
  \code{\link{lda_cgs_perplexity}},
  \code{\link{lda_fgs_BF_perplexity}},
  \code{\link{lda_fgs_perplexity}},
  \code{\link{lda_fgs_st_perplexity}}
}
