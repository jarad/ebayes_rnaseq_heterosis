data {
  int<lower=1> S;
  int<lower=0> count[S]; 
  matrix[S,3] X;
  real eta_phi;
  real eta_alpha;
  real eta_delta;
  real eta_psi;
  real sigma_phi;
  real sigma_alpha;
  real sigma_delta;
  real sigma_psi;
  vector[S] c;                     // lane sequencing depth
}
parameters {
  real phi;
  real alpha;
  real delta;
  real psi;          
}
transformed parameters {
  vector[3] pad;

  pad[1] <- phi;
  pad[2] <- alpha;
  pad[3] <- delta;
}
model {
  phi   ~ normal(            eta_phi,   sigma_phi);
  alpha ~ double_exponential(eta_alpha, sigma_alpha); // Laplace
  delta ~ double_exponential(eta_delta, sigma_delta); // Laplace
  psi   ~ normal(            eta_psi,   sigma_psi);

  count ~ neg_binomial_2_log(X*pad+c, 1/exp(psi));
}
