data {
  int<lower=1> S;
  int<lower=0> count[S]; 
  int<lower=1,upper=2> variety[S];
  real eta_phi;
  real eta_alpha;
  real eta_psi;
  real sigma_phi;
  real sigma_alpha;
  real sigma_psi;
  vector[S] c;                     // lane sequencing depth
}
transformed data {
  matrix[S,2] X;
  for (s in 1:S) {
    if (variety[s] == 1) { X[s,1] <-  1; X[s,2] <- -1; }
    if (variety[s] == 2) { X[s,1] <-  1; X[s,2] <-  1; }
  }
}
parameters {
  real phi;
  real alpha;
  real psi;          
}
transformed parameters {
  vector[2] pa;

  pa[1] <- phi;
  pa[2] <- alpha;
}
model {
  phi   ~ normal(            eta_phi,   sigma_phi);
  alpha ~ double_exponential(eta_alpha, sigma_alpha); // Laplace
  psi   ~ normal(            eta_psi,   sigma_psi);

  count ~ neg_binomial_2_log(X*pa+c, 1/exp(psi));
}
