data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=1,upper=G> gene[N];
  real y[N];
  real mu;
  real<lower=0> tau;
  real<lower=0> sigma;
}

parameters {
  real theta[G];
}

model {
  for (g in 1:G) theta[g] ~ normal(mu, tau);
  for (n in 1:N) y[n] ~ normal(theta[gene[n]], sigma);
}
