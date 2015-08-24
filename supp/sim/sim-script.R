library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

source("args.R")

source("../common/get_hyperparameters.R")
source("../common/single_gene_analysis.R")

# Compile stan model
if (m == 'laplace') model = stan_model("../common/laplace.stan")
if (m == 'normal' ) model = stan_model("../common/normal.stan")


data_file = paste0("data/sim-", r, "-", i,'.rds')

d = as.data.frame(readRDS(data_file)); n = ncol(d)/3
names(d) = c(paste("B73_",1:n, sep=''), paste('Mo17_',1:n, sep=''), paste('B73xMo17_',1:n, sep=''))
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# trim genes for average count greater than 1
# d = d[which(rowMeans(d)>1),] 
hyperparameters = get_hyperparameters(d, variety, method=m)



# ######################################
# Run individual gene analyses
######################################

# This will take a long time, so parallelize if possible
if (parallel <- require(doMC)) {
  registerDoMC()
} else if (parallel <- require(doParallel)) {
  cl = makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
}

# Divert Stan out put to temporary file
tmp_file = paste0("Rout/", m, "/sim-", r, "-", i,'.tmp')
sink(file=tmp_file)

# Run MCMC analysis on each gene
analysis = adply(d,
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene', 
                 .parallel = parallel,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))

# Remove temporary file
sink()
unlink(tmp_file)


rownames(analysis) = rownames(d)


results = analysis[,c("prob_LPH","prob_HPH","effectiveSize")]

saveRDS(results, file=paste0("results/", m, "/sim-", r, "-", i,'.rds'))


q(ifelse(interactive(), "ask","no"))
