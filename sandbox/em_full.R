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

G = 1000
d = sim_heterosis_data(G=G)


# Set initial values for hyperparameters
gene_names = c("phi","alpha","delta","psi")
hyper = with(d$hyperparameters, as.list(c(location, scale)))
names(hyper) = c(paste("eta_",   gene_names, sep=''),
                 paste("sigma_", gene_names, sep=''))
n_hyper = 8

truth = data.frame(variable=names(hyper), value=with(d$hyperparameter, c(location,scale)))

# 
n_iter = 100
hyper_keep = matrix(NA, nrow=n_iter, ncol=n_hyper)

cat("MCEM with", G, "genes.\n", file="em_full.Ro")
cat("iteration", names(hyper), "iteration_time\n", sep="")
# Run EM
for (i in 1:n_iter) {
  start_time = proc.time()
  
  # MCMC for gene specific parameters
  mcmc = dlply(d$data, .(gene), function(x) {
    sampling(m, c(list(S=nrow(x), variety=x$parent, count=x$y, c=rep(1, nrow(x))), hyper), gene_names, iter=100+(i+2)^2, chains=1)
  }, .parallel=TRUE)
  for(j in 1:length(mcmc)) attr(mcmc[[j]],"name") <- names(mcmc)[j]

  # Get initial values (for next step)
  
  
  # Extract samples to calculate M step
  samps = ldply(mcmc, function(x) {
    tmp = extract(x, gene_names)
    data.frame(gene = attr(x, 'name'), phi = tmp$phi, alpha = tmp$alpha, delta = tmp$delta, psi = tmp$psi)
  })
  
  # Calculate summary statistics
  sm = ddply(samps, .(gene), summarize,
             phi_mean = mean(phi),
             phi_var  = var(phi),
             psi_mean = mean(psi),
             psi_var  = var(psi),
             alpha_median = median(alpha),
             delta_median = median(delta), 
             .parallel=TRUE)

  # Update hyperparameters
  hyper$eta_phi = mean(sm$phi_mean)
  hyper$eta_psi = mean(sm$psi_mean)
  hyper$sigma_phi = sqrt(mean(sm$phi_var + (sm$phi_mean-hyper$eta_phi)^2))
  hyper$sigma_psi = sqrt(mean(sm$psi_var + (sm$psi_mean-hyper$eta_psi)^2))
  
  hyper$eta_alpha = median(sm$alpha_median)
  hyper$eta_delta = median(sm$delta_median)
  sm2 = ddply(samps, .(gene), summarize,
              alpha_exabsd = mean(abs(alpha-hyper$eta_alpha)),
              delta_exabsd = mean(abs(delta-hyper$eta_delta)), .parallel=TRUE)
  hyper$sigma_alpha = mean(sm2$alpha_exabsd)
  hyper$sigma_delta = mean(sm2$delta_exabsd)
  
  # Save hyperparameters
  for (j in 1:n_hyper) hyper_keep[i,j] = hyper[[j]]
  
  save(hyper_keep, file="em_full.RData")
  
  cat(i, unlist(hyper), proc.time()[3]-start_time[3], "\n",
      file="em_full.Ro", append=TRUE)
}

# Plot hyperparameters
em = data.frame(iteration=1:n_iter, hyper_keep)
names(em)[-1] = names(hyper)

ggplot(melt(em, id.var='iteration'), aes(x=iteration, y=value)) + 
  geom_point() +
  facet_wrap(~variable, scale='free', nrow=2) +
  geom_hline(data=truth, aes(yintercept=value), col='red')


