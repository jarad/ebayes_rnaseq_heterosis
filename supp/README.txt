data.csv:   the B73-Mo17 data set
model.stan: the stan model file

script.R:     an R script to run the data analysis
figs.R: an R script to create figures (depends on script.R)

get_hyperparameters.R: function to obtain hyperparameters
single_gene_analysis.R: function to run MCMC for a gene conditional on hyperparameters
