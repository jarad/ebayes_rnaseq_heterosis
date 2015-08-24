library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("../common/get_hyperparameters.R")
source("../common/single_gene_analysis.R")

##################################################
# Get hyperparameters
##################################################
d = read.csv("data/data.csv", row.names=1)
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# trim genes for average count greater than 1
#d = d[which(rowMeans(d)>1),] 
hyperparameters = get_hyperparameters(d, variety, method='laplace')



##################################################
# Perform individual gene analsysis (in parallel)
##################################################

# Compile stan model
model = stan_model("../common/laplace.stan")

# This will take a long time, so parallelize if possible
require(doMC)
registerDoMC(5)

sink("script.temp")

analysis = adply(d,
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene',
                 .parallel = TRUE,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))

sink()
unlink("script.temp")

rownames(analysis) = rownames(d)
saveRDS(analysis, file="results/laplace-parallel.rds")

q(ifelse(interactive(), "ask", "no"))
