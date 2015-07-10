data {
  int<lower=1> S;
  int<lower=0> count[S]; 
  int<lower=1,upper=3> variety[S];
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
transformed data {
  matrix[S,3] X;
  for (s in 1:S) {
    if (variety[s] == 1) { X[s,1] <-  1; X[s,2] <- -1; X[s,3] <-  0; }
    if (variety[s] == 2) { X[s,1] <-  1; X[s,2] <-  1; X[s,3] <-  0; }
    if (variety[s] == 3) { X[s,1] <-  1; X[s,2] <-  0; X[s,3] <-  1; }
  }
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
generated quantities {
  int<lower=0, upper=1> LPH;
  int<lower=0, upper=1> HPH;
  
  LPH <-      delta  < -fabs(alpha);
  HPH <- fabs(delta) >  fabs(alpha);
}
