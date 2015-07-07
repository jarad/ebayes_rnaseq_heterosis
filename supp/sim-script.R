library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("get_hyperparameters.R")
source("single_gene_analysis.R")

# Compile stan model
model = stan_model("model.stan")

rds_files = system("ls -S -r ~/temp/rna-seq-stan/simulations/sim1/datasets/*.rds", TRUE)
results = list()
for (rds in rds_files) {
  
  d = as.data.frame(readRDS(rds)); n = ncol(d)/3
  names(d) = c(paste("B73_",1:n, sep=''), paste('Mo17_',1:n, sep=''), paste('B73xMo17_',1:n, sep=''))
  variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))
  
  # trim genes for average count greater than 1
  # d = d[which(rowMeans(d)>1),] 
  hyperparameters = get_hyperparameters(d, variety)

  
  
  ######################################
  # Run individual gene analyses
  ######################################
  
  # This will take a long time, so parallelize if possible
  if (parallel <- require(doMC)) {
    registerDoMC()
  } else if (parallel <- require(doParallel)) {
    cl = makeCluster(detectCores(), type='PSOCK')
    registerDoParallel(cl)
  }
  
  analysis = adply(d,
                   1,
                   function(x) single_gene_analysis(x),
                   .id = 'gene', 
                   .parallel = parallel,
                   .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))
  
  rownames(analysis) = rownames(d)
  

  results[[rds]] = analysis[,c("prob_LPH","prob_HPH","effectiveSize")]
  saveRDS(results, file='sim-results.RDS')
}

q(ifelse(interactive(), "ask","no"))
