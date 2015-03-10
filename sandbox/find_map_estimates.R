library(Rcpp)
library(reshape2)

source("sim_heterosis_data.R")
source("Laplace.R")

# Get data
all = sim_heterosis_data(G=100)

# Reshape data (put this in a function?)
count = dcast(all$data[,c("gene","sample","y")], gene ~ sample, value.var='y')
variety = all$data$parent[1:ncol(count)-1] # hack


hyperparameters = melt(all$hyperparameters)

sourceCpp("monte_carlo_integral.cpp")

log_like = function(pi) {
  -monte_carlo_integral(t(as.matrix(count[,-1])), variety-1, pi[1:4], exp(pi[5:8]), 1000)
}

o = optim(rep(0,8), log_like)
