library(Rcpp)
library(reshape2)

source("sim_heterosis_data.R")
source("Laplace.R")

# Get data
all = sim_heterosis_data(G=3)

# Reshape data (put this in a function?)
count = dcast(all$data[,c("gene","sample","y")], gene ~ sample, value.var='y')
variety = all$data$parent[1:ncol(count)-1] # hack


hyperparameters = melt(all$hyperparameters)

sourceCpp("monte_carlo_integral.cpp")

integral = rep(NA,20)
for (i in 1:length(integral))
  integral[i] = monte_carlo_integral(as.matrix(count[,-1]), variety-1, 
                     with(hyperparameters, value[variable=="location"]), 
                     with(hyperparameters, value[variable=="scale"]), 
                     4)

