data {
  int<lower=1> S;
  int<lower=0> count[S]; 
  int<lower=1,upper=3> variety[S];
  real eta_phi;
  real eta_alpha;
  real eta_delta;
  real sigma_phi;
  real sigma_alpha;
  real sigma_delta;
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
