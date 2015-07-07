library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("get_hyperparameters.R")
source("single_gene_analysis.R")

##################################################
# Get hyperparameters
##################################################
d = read.csv("data.csv", row.names=1)
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# trim genes for average count greater than 1
#d = d[which(rowMeans(d)>1),] 

hyperparameters = get_hyperparameters(d, variety)



##################################################
# Perform individual gene analsysis (in parallel)
##################################################

# Compile stan model
model = stan_model("model.stan")

# This will take a long time, so parallelize if possible
if (parallel <- require(doMC)) {
  registerDoMC()
} else if (parallel <- require(doParallel)) {
  cl = makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
}

wh = 1:200
analysis = adply(d[wh,],
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene',
                 .parallel = parallel,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))

