library(rstan)
library(plyr)

# Compile models
models = list()
models$de_poisson  = stan_model('stan_models/de_poisson.stan')
models$de_negbin   = stan_model('stan_models/de_negbin.stan')
models$het_poisson = stan_model('stan_models/het_poisson.stan')
models$het_negbin  = stan_model('stan_models/het_negbin.stan')


# calculates hyper parameters for normal distribution
normal_hyper = function(samples, learn_mean) {
  eta = ifelse(learn_mean, mean(samples, 0))
  sigma = sqrt(mean((samples-eta)^2))
  return(list(eta=eta,sigma=sigma))
}

# calculates hyper parameters for laplace distribution
laplace_hyper = function(samples, learn_mean) {
  eta = ifelse(learn_mean, median(samples,0))
  sigma = mean(abs(samples-eta))
  return(list(eta=eta,sigma=sigma))
}

# MCEM algorithm for 
em_de_poisson = function(dat, niter, learn_alpha_mean = FALSE, alpha_dist = 'normal') {
  parameter_names = c("phi","alpha")
  
  # Build X matrix
  one_gene = dat[dat$gene==1,]
  X = cbind(1, 2*(one_gene$variety-1.5))
  
  # Initial values
  initial = ddply(dat, .(gene), function(x) {
    m = glm(count~0+X, data=x, family="poisson")
    data.frame(phi = coef(m)[1], alpha = coef(m)[2])
  })
  hyper = list()
  hyper$iteration = 0
  hyper$eta_phi = mean(initial$phi)
  hyper$sigma_phi = sd(initial$phi)
  if (alpha_dist == 'normal') {
    hyper$eta_alpha = ifelse(learn_eta_alpha, mean(initial$alpha), 0)
    hyper$sigma_alpha = sqrt(mean((initial$alpha-hyper$eta_alpha)^2))
  } else if (alpha_dist == 'laplace') {
    hyper$eta_alpha = ifelse(learn_eta_alpha, median(initial$alpha), 0)
    hyper$sigma_alpha = mean(abs(initial$alpha-hyper$eta_alpha))
  } else {
    stop(paste('Alpha distribution',alpha_dist,'not implemented.'))
  }
  
  # Save structure
  hyper_keep = as.data.frame(hyper)
  
  # MCEM
  for (i in 1:niter) {
    hyper$iteration = i
    # (Monte Carlo) Expectation
    mcmc = dlply(dat, .(gene), function(x) {
      sampling(models$de_poisson, c(list(S=nrow(x), 
                         variety = x$variety, 
                         count=x$count, 
                         c=rep(0, nrow(x))), # fix this
                    hyper), 
               parameter_names, 
               iter=10^2+(i+2)^2, 
               chains=4)
    }, .parallel=TRUE)
    for(j in 1:length(mcmc)) attr(mcmc[[j]],"name") <- names(mcmc)[j]
    
    # Maximization
    # Obtain MCMC samples
    samps = ldply(mcmc, function(x) {
      tmp = extract(x, parameter_names)
      data.frame(gene  = attr(x, 'name'), 
                 phi   = tmp$phi, 
                 alpha = tmp$alpha)
    }, .parallel=TRUE)
    
    
    # Calculate gene specific summary statistics
    sm = ddply(samps, .(gene), summarize,
               phi_mean = mean(phi),
               phi_var  = var(phi),
               alpha_median = median(alpha),
               .parallel=TRUE)
    
    # phi
    tmp = normal_hyper(samps$phi, TRUE)
    hyper$eta_phi = tmp$eta
    hyper$sigma_phi = tmp$sigma
    
    # alpha
    if (alpha_dist=='normal' ) {
      tmp = normal_hyper(samps$alpha, learn_alpha_mean)
    } else if (alpha_dist=='laplace') {
      tmp = laplace_hyper(samps$alpha, learn_alpha_mean)
    }
    hyper$eta_alpha = tmp$eta
    hyper$sigma_alpha = tmp$sigma
    
    # Save
    hyper_keep = rbind(hyper_keep, as.data.frame(hyper))
  }
    
  # Return
  hyper_keep
}
