ldamcmc
=======

This **R** package implements several Markov chain Monte Carlo (MCMC) algorithms for the latent Dirichlet allocation (LDA) model. 


Dependencies (*R* packages) 
----------------------------

* **Rcpp**
* **RcppArmadillo** that uses the **Armadillo** C++ package 
* **MCMCpack** 
* **lattice**

Installation Guide 
------------------

* Download the package source from [Git Download Link](https://github.com/clintpgeorge/ldamcmc/archive/master.zip)
* Unzip the dowloaded file and rename the folder **ldamcmc-master** to **ldamcmc** 
* To install this package, run ***R CMD INSTALL ldamcmc*** on the commandline 
* To uninstall this package, run ***R CMD REMOVE ldamcmc*** on the commandline 

References
----------

* George, C.P. and Doss, H. (2015). Principled Selection of Hyperparameters in the Latent Dirichlet Allocation Model. Technical Report, Department of Computer and Information Science and Engineering, University of Florida
* Chen, Z. and Doss, H. (2015). Inference for the Number of Topics in the Latent Dirichlet Allocation Model via Bayesian Mixture Modelling. Technical Report, Department of Statistics, University of Florida.
