ldamcmc
=======

This **R** package implements several Markov chain Monte Carlo (MCMC) algorithms for the latent Dirichlet allocation (LDA) model. 


Dependencies
----------------------------

This package uses the following R packages, which are already included in this R package.   
* **Rcpp**
* **RcppArmadillo** based on the **Armadillo** C++ package 
* **MCMCpack** 
* **lattice**

Installation Guide 
------------------

* Download the package source from [Git Download Link](https://github.com/clintpgeorge/ldamcmc/archive/master.zip)
* Unzip the dowloaded file and rename the folder **ldamcmc-master** to **ldamcmc** 
* To install this package, run ```R CMD INSTALL ldamcmc``` on the commandline 
* To uninstall this package, run ```R CMD REMOVE ldamcmc``` on the commandline 

References
----------

1. George, C.P. and Doss, H. (2015). Principled Selection of Hyperparameters in the Latent Dirichlet Allocation Model. Technical Report, Department of Computer and Information Science and Engineering, University of Florida
2. Chen, Z. and Doss, H. (2015). Inference for the Number of Topics in the Latent Dirichlet Allocation Model via Bayesian Mixture Modelling. Technical Report, Department of Statistics, University of Florida.

Acknowledgements
----------------

This research is partially supported by generous contributions from the International Center for Automated Research at the University of Florida Levin College of Law. 

The author would like to express sincere thanks to Dr. Hani Doss for all the insightful comments and advice on this research. The author also would like to thank Dr. Joseph N. Wilson, Zhe Chen, and Wei Xia for the valuable discussions and critics that helped throughout the development of this project.

