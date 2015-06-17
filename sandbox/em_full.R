library(rstan)
library(reshape2)
library(doMC)
#library(dplyr)
library(ggplot2)
library(plyr)

source("Laplace.R")
source("sim_heterosis_data.R")

registerDoMC(cores=8)

m = stan_model("sg_model_pad_dexp.txt")

set.seed(20150408)
G = 100
d = sim_heterosis_data(G=G)

# Construct design
X = cbind(1, rep(c(-1,1,0), each=4), rep(c(0,0,1), each=4))
#X = cbind(1, rep(c(0,1,0), each=4), rep(c(0,0,1), each=4))

# Use edgeR
source("edgeR_est.R")
d_w = dcast(d$data[,c("gene","parent","sample","y")], gene~parent+sample, value.var='y')
edgeR = edgeR_est(d_w[,-1], substr(names(d_w)[-1],1,1), rownames(d_w))

# Set initial values for hyperparameters
#gene_names = c("phi","alpha","delta","psi")
gene_names = c("phi","alpha","delta")
hyper = with(edgeR$hyperparameters[1:3,], as.list(c(location, scale)))
names(hyper) = c(paste("eta_",   gene_names, sep=''),
                 paste("sigma_", gene_names, sep=''))
hyper$eta_alpha = 0 
hyper$eta_delta = 0
n_hyper = length(hyper)

truth = data.frame(variable=names(hyper), value=with(d$hyperparameter, c(location,scale)))

# 
n_iter = 10
hyper_keep = matrix(NA, nrow=n_iter, ncol=n_hyper)

if (!interactive()) {
  cat("MCEM with", G, "genes.\n", file="em_full.Ro")
  cat("iteration", names(hyper), "iteration_time\n", sep=", ", file="em_full.Ro", add=TRUE)
}



# Run EM
for (i in 1:n_iter) {
  start_time = proc.time()
  
  # MCMC for gene specific parameters
  mcmc = dlply(d$data, .(gene), function(x) {
    sampling(m, c(list(S=nrow(x), 
                       X = X, 
                       count=x$y, 
                       c=rep(0, nrow(x))), 
                  hyper), 
             c(gene_names), 
             iter=10^2+(i+2)^2, 
             chains=4)
  }, .parallel=TRUE)
  for(j in 1:length(mcmc)) attr(mcmc[[j]],"name") <- names(mcmc)[j]

  # Get initial values (for next iteration)
  
  
  # Extract samples to calculate M step
  samps = ldply(mcmc, function(x) {
    tmp = extract(x, gene_names)
    data.frame(gene  = attr(x, 'name'), 
               phi   = tmp$phi, 
               alpha = tmp$alpha, 
               delta = tmp$delta)
#               psi   = tmp$psi
  })
  
  # Calculate summary statistics
  sm = ddply(samps, .(gene), summarize,
             phi_mean = mean(phi),
             phi_var  = var(phi),
#             psi_mean = mean(psi),
#             psi_var  = var(psi),
             alpha_median = median(alpha),
             delta_median = median(delta),
             .parallel=TRUE)
  
  # Update hyperparameters
  hyper$eta_phi = mean(samps$phi)
#hyper$eta_alpha = mean(samps$alpha)
#hyper$eta_delta = mean(samps$delta)
#  hyper$eta_psi = mean(samps$psi)
  hyper$sigma_phi = sqrt(mean((samps$phi-hyper$eta_phi)^2))
#  hyper$sigma_psi = sqrt(mean((samps$psi-hyper$eta_psi)^2))
  hyper$sigma_alpha = mean(abs(samps$alpha))
  hyper$sigma_delta = mean(abs(samps$delta))
   
#   hyper$eta_alpha = median(sm$alpha_median) 
#   hyper$eta_delta = median(sm$delta_median) 
#   sm2 = ddply(samps, .(gene), summarize,
#               alpha_exabsd = mean(abs(alpha-hyper$eta_alpha)),
#               delta_exabsd = mean(abs(delta-hyper$eta_delta)), .parallel=TRUE)
#   hyper$sigma_alpha = mean(sm2$alpha_exabsd)
#   hyper$sigma_delta = mean(sm2$delta_exabsd)
  

  
  
  # Save hyperparameters
  hyper_keep[i,] = unlist(hyper)
  
  if (!interactive()) {
    save(hyper_keep, file="em_full.RData")
  
    cat(i, unlist(hyper), proc.time()[3]-start_time[3], "\n",
        file="em_full.Ro", append=TRUE)
  } else {
    cat(i, unlist(hyper),"\n")
  }
}

# Plot hyperparameters
em = data.frame(iteration=1:n_iter, hyper_keep)
names(em)[-1] = names(hyper)

ggplot(melt(em, id.var='iteration'), aes(x=iteration, y=value)) + 
  geom_point() +
  facet_wrap(~variable, scale='free', nrow=2) 
+
  geom_hline(data=truth[1:n_hyper,], aes(yintercept=value), col='red')


