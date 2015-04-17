data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=1,upper=G> gene[N];
  int<lower=0> y[N];
  real mu;
  real<lower=0> tau;
  real<lower=0> sigma;
}

parameters {
  real theta[G];
}

model {
  theta ~ normal(mu, tau);
  for (n in 1:N) y[n] ~ neg_binomial_2(theta[gene[n]], sigma);
}
