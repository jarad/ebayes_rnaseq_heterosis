library(reshape2)
library(ggplot2)
library(plyr)
library(lme4)
library(distr)
library(rstan)

# Compile Stan models
normal_em = stan_model("normal_em.stan")
laplace_em = stan_model("laplace_em.stan")
nb_normal = stan_model("nb_normal.stan")


####################################################
# Normal random effects example
####################################################
set.seed(1)

G = 1000
n = 10

truth = data.frame(mu=0, tau=1, sigma=1)

d = ddply(data.frame(gene=1:G), .(gene), function(x) {
  theta = rnorm(1, truth$mu, truth$tau)
  data.frame(theta = theta, y = rnorm(n, theta, truth$sigma))
})

sm = ddply(d, .(gene), summarize, mean=mean(y), sd=sd(y))

# Truth
o = lmer(y~(1|gene), d, REML=FALSE)
(truth = data.frame(parameter=c("mu","tau","sigma"), value=c(fixef(o), as.data.frame(VarCorr(o))$sdcor)))
# alternatively
mean(d$y)           # == truth[1]

sqrt(mean(sm$sd^2)) # == truth[3]

# EM
mu = 0
tau = 1
sigma = truth$value[truth$parameter=='sigma'] # only writing the EM for mu and tau

n_iter = 10
mu_em = rep(mu,10)
tau_em= rep(tau,10)

for (i in 2:n_iter) {
  theta_var  = 1/(n/sigma^2+1/tau^2)
  theta_mean = theta_var*(n*sm$mean/sigma^2+mu/tau^2)
    
  # Maximization
  mu = mean(theta_mean)
  tau = sqrt(mean(theta_var + (theta_mean-mu)^2)) # add the variance and deviation from the mean
  
  # Save values
  mu_em[i] = mu
  tau_em[i] = tau
}

em = data.frame(iteration=1:length(mu_em),
                mu = mu_em,
                tau = tau_em,
                method = "em")

#######################################
# MCEM
#######################################
# Initial values
mu =0 
tau = 1

n_iter = 10
mu_em_stan  = rep(mu,10)
tau_em_stan = rep(tau,10)

pre_dat = list(N=nrow(d), G=G, y=d$y, gene=d$gene, sigma=truth$value[truth$parameter=='sigma'] )

for (i in 2:n_iter) {
  # Run MCMC
  r = sampling(normal_em, c(pre_dat, list(mu=mu, tau=tau)))
  theta = extract(r, "theta")$theta
  theta_mean = colMeans(theta) # aaply(theta, 2, mean)
  theta_var  = aaply(theta, 2, var)
  
  # Maximization
  mu = mean(theta_mean)
  tau = sqrt(mean(theta_var + (theta_mean-mu)^2)) # add the variance and deviation from the mean
  
  # Save values
  mu_em_stan[i] = mu
  tau_em_stan[i] = tau
}


mc_em = data.frame(iteration=1:length(mu_em_stan),
                mu = mu_em_stan,
                tau = tau_em_stan,
                method = "mc_em")

ggplot(melt(rbind(em,mc_em), id=c('iteration','method'), variable.name='parameter'), aes(x=iteration, y=value, color=method)) + 
  geom_point() + 
  geom_hline(data=truth[1:2,], aes(yintercept=value)) + 
  facet_wrap(~parameter, scales='free') 

##############################################
# Laplace random effects example
##############################################
set.seed(1)

G = 1000
n = 9

truth = list(mu=0, tau=1, sigma=1)

d = ddply(data.frame(gene=1:G), .(gene), function(x) {
  theta = truth$mu + r(DExp(rate=truth$tau))(1)
  data.frame(theta = theta, y = rnorm(n, theta, truth$sigma))
})

# sm = ddply(d, .(gene), summarize, 
#            mean   = mean(y)
#            median = median(y), 
#            mad=mad(y, constant=1), # Our random effects are not normal
#            ) 

# Truth
# ?? stan ??

# Initial values
mu = 0 
tau = 1

n_iter = 10
mu_em_stan  = rep(mu,10)
tau_em_stan = rep(tau,10)

pre_dat = list(N=nrow(d), G=G, y=d$y, gene=d$gene, sigma=1)

for (i in 2:n_iter) {
  # Run MCMC
  r = sampling(laplace_em, c(pre_dat, list(mu=mu, tau=tau)))
  theta = extract(r, "theta")$theta
  theta_median = aaply(theta, 2, median) # aaply(theta, 2, mean)
  
  # Maximization
  mu = median(theta_median)
  
  # Calculate gene expected absolute deviations
  theta_edev = aaply(theta,2, function(x) mean(abs(x-mu)))
  tau = sqrt(mean(theta_edev)) 
  
  # Save values
  mu_em_stan[i] = mu
  tau_em_stan[i] = tau
}


mc_em = data.frame(iteration=1:length(mu_em_stan),
                   mu = mu_em_stan,
                   tau = tau_em_stan,
                   method = "mc_em")

ggplot(melt(mc_em_laplace, id=c('iteration','method'), variable.name='parameter'), aes(x=iteration, y=value, color=method)) + 
  geom_point() + 
#  geom_hline(data=truth[1:2,], aes(yintercept=value)) + 
  facet_wrap(~parameter, scales='free') 




###################################
# Negative binomial model
###################################

# Normal prior
###################################
set.seed(1)

G = 1000
n = 10

truth = data.frame(mu=0, tau=1, sigma=1)

d = ddply(data.frame(gene=1:G), .(gene), function(x) {
  theta = rnorm(1, truth$mu, truth$tau)
  data.frame(theta = theta, y = rnbinom(n, mu=exp(theta), size=truth$sigma))
})

# Truth
o = glmer.nb(y~(1|gene), data=d, nAGQ=1)
(truth = data.frame(parameter=c("mu","tau","sigma"), value=c(fixef(o), as.data.frame(VarCorr(o))$sdcor)))


# MCEM
#######################################
# Initial values
mu =0 
tau = 1

n_iter = 10
nb_mu_em  = rep(mu,10)
nb_tau_em = rep(tau,10)

pre_dat = list(N=nrow(d), G=G, y=d$y, gene=d$gene, sigma=truth$value[truth$parameter=='sigma'] )

for (i in 2:n_iter) {
  # Run MCMC
  r = sampling(normal_em, c(pre_dat, list(mu=mu, tau=tau)))
  theta = extract(r, "theta")$theta
  theta_mean = colMeans(theta) # aaply(theta, 2, mean)
  theta_var  = aaply(theta, 2, var)
  
  # Maximization
  mu = mean(theta_mean)
  tau = sqrt(mean(theta_var + (theta_mean-mu)^2)) # add the variance and deviation from the mean
  
  # Save values
  nb_mu_em[i] = mu
  nb_tau_em[i] = tau
}


mc_em = data.frame(iteration=1:length(mu_em_stan),
                   mu = mu_em_stan,
                   tau = tau_em_stan,
                   method = "mc_em")

ggplot(melt(mc_em, id=c('iteration','method'), variable.name='parameter'), aes(x=iteration, y=value, color=method)) + 
  geom_point() + 
  geom_hline(data=truth[1:2,], aes(yintercept=value)) + 
  facet_wrap(~parameter, scales='free') 

