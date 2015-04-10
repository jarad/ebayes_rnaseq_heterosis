data {
  int<lower=1> S;
  int<lower=0> count[S]; 
  matrix[S,3] X;
  real eta_phi;
  real eta_alpha;
  real eta_delta;
  real sigma_phi;
  real sigma_alpha;
  real sigma_delta;
  vector[S] c;                     // lane sequencing depth
}
parameters {
  real phi;
  real alpha;
  real delta;        
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

  count ~ poisson_log(X*pad + c);
}
