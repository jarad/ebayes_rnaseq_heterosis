library(plyr)

G = 10
v = 3
nv = rep(4,3)

rlaplace = function(n, mean, scale) {
  u = runif(n, -.5, .5)
  mean - scale * sign(u) * log(1-2*abs(u))
}

# Priors
prior = list(eta_phi     = 1,
             eta_alpha   = 0,
             eta_delta   = 0,
             eta_psi     = -2,
             sigma_phi   = 1,
             sigma_alpha =0.1,
             sigma_delta = 0.01,
             sigma_psi   = 0.1)

# Parameters
parameters = with(prior, list(phi   = rnorm   (G, eta_phi,   sigma_phi  ),
                              alpha = rlaplace(G, eta_alpha, sigma_alpha),
                              delta = rlaplace(G, eta_delta, sigma_delta),
                              psi   = rnorm   (G, eta_psi,   sigma_psi  )))

# Data
d = ddply(data.frame(gene=1:G), .(gene), function(x) {
  mutate(data.frame(parent = rep(1:3, times = nv)),
         sample = 1:sum(nv),
         mu = with(parameters, c(phi[g]+alpha[g], phi[g]-alpha[g], phi[g]+delta[g])),
         y = # need negative binomial with mu parameterization
  )
  
})
