Markov Chain Monte Carlo Algorithms for the Latent Dirichlet Allocation Model
=============================================================================

This **R** package, **ldamcmc**, implements several Markov chain Monte Carlo (MCMC) algorithms for the latent Dirichlet allocation (LDA) model. This includes: 

* The augmented collapsed Gibbs sampling (ACGS, Griffiths and Steyvers 2004, George and Doss 2015) algorithm
* The full Gibbs sampling (FGS, George and Doss 2015) algorithm
* The serial tempering (George and Doss 2015, Geyer 2011) algorithm 
* Hyperparameter selection in the LDA model (George and Doss 2015) 
* Posterior predictive checking (PPC, Chen and Doss 2015)

For package documentation run 

``` help("ldamcmc") ```

in an R console. All major functions and datasets are documented and linked to the package index. Raw data files for each dataset are available in the **data-raw** folder. To load raw data see ``` demo/load_raw_data.R  ```.    

To see all demo R scripts available in this package, run 

``` demo(package="ldamcmc") ```

in an R console. Some scripts can be executed via running  

``` demo(file-name, package="ldamcmc") ```

in an R console. The rest of them may require commandline arguments for execution. Please see the documentation provided in each script before execution.    

Authors
----------------------------
* [Clint P. George](http://www.cise.ufl.edu/~cgeorge) (Please contact for questions and comments)
* [Hani Doss](http://www.stat.ufl.edu/~doss) 

Dependencies
----------------------------

This package uses the following R packages, which are already included in this R package.   
* **Rcpp**
* **RcppArmadillo** based on the **Armadillo** C++ package 
* **MCMCpack**, which may require packages such as **coda** and **MASS** 
* **lattice**

Installation Guide 
------------------

* Download the package source from [Git Download Link](https://github.com/clintpgeorge/ldamcmc/archive/master.zip)
* Unzip the dowloaded file and rename the folder **ldamcmc-master** to **ldamcmc** 
* To install **ldamcmc** run ```R CMD INSTALL ldamcmc``` on the commandline 
* To uninstall **ldamcmc** run ```R CMD REMOVE ldamcmc``` on the commandline 

References
----------

1. Blei, D. M., Ng, A. Y. and Jordan, M. I. (2003). Latent Dirichlet 
allocation. Journal of Machine Learning Research 3 993-1022.
2. Chen, Z. and Doss, H. (2015). Inference for the number of topics in the 
latent Dirichlet allocation model via Bayesian mixture modelling. Tech. rep., 
Department of Statistics, University of Florida.
3. George, C.P. and Doss, H. (2015). Principled Selection of Hyperparameters 
in the Latent Dirichlet Allocation Model. Tech. rep., Department of 
Computer and Information Science and Engineering, University of Florida 
4. Geyer, C. J. (2011). Importance sampling, simulated tempering, and 
umbrella sampling. In Handbook of Markov Chain Monte Carlo (S. P. Brooks, A. 
E. Gelman, G. L. Jones and X. L. Meng, eds.). Chapman & Hall/CRC, Boca Raton, 
295-311.
5. Griffiths, T. L. and Steyvers, M. (2004). Finding scientific topics. 
Proceedings of the National Academy of Sciences 101 5228-5235.

Acknowledgements
----------------

This research is partially supported by generous contributions from the [International Center for Automated Research](http://www.law.ufl.edu/academics/institutes/icair) (ICAIR) at the University of Florida Levin College of Law. 

The authors would like to thank Dr. Joseph N. Wilson, Zhe Chen, and Wei Xia for the valuable discussions and critics that helped throughout the development of this project.

