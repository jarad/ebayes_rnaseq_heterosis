common/ a directory for common files
  empirical_estimate_function.R: a function to implement Ji et al.
  get_hyperparameters.R: function to obtain hyperparameters
  normal.stan:  stan model file with normal distributions
  laplace.stan: stan model file with Laplace distributions
  single_gene_analysis.R: function to run MCMC for a gene conditional on hyperparameters

eric: a directory for Eric's adjusted work (don't need this in the supplement)

real/: a directory for the real data analysis
  data/: a directory for the data
    data.csv:   the B73-Mo17 data set
  figs.R: an R script to create figures (depends on script.R)
  ji.R: an R script to run Ji et al. 
  script.R:     an R script to run the data analysis


sim/: a directory simulation study
  data/: a directory for simulated data
  results/: a directory for results
    ji/: Ji et al. results
    normal/: eBayes with normal distribution results
    laplace/: eBayes with Laplace distribution results
  args.R: a script to take care of command line arguments
  create_Makefile.R: an R script to create the Makefile
  ji-sim-script.R: a script to perform Ji et al. on simulated data
  sim-script.R: a script to perform eBayes analysis of simulation data

