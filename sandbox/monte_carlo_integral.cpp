#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

double logsumexp(double x[], int n) {
  double max_val = x[0], sum = 0.0;
  int i;
  
  for (i = 1; i< n; i++)
    if (x[i] > max_val) 
      max_val = x[i];
      
  for (i=0; i<n; i++)
    sum += exp(x[i] - max_val);
    
  return log(sum) + max_val;
}

// [[Rcpp::export]]
double monte_carlo_integral(IntegerMatrix count, 
                            IntegerVector variety, 
                            NumericVector hyperparameters,
                            int n_sims) {
  // count is the observed count (gene x sample)
  // variety is the variety identifier for the columns
  // hyperparameters are the values for the hyperparameters
  // sims number of Monte Carlo simulations to use
  int G = count.nrow();
  int n = count.ncol();
  
  // argument error handling
  if (variety.size() != n) throw std::range_error("Length of variety and number of columns of count do not agree.");
  if (hyperparameters.size() != 8) throw std::range_error("Hyperparameters must have length 8.");
  
  // set up hyperparameters
  double phi_location   = hyperparameters[0],
         alpha_location = hyperparameters[1],
         delta_location = hyperparameters[2],
         psi_location   = hyperparameters[3],
         phi_scale      = hyperparameters[4],
         alpha_scale    = hyperparameters[5], 
         delta_scale    = hyperparameters[6],
         psi_scale      = hyperparameters[7];
         
  double integral = 0, partial_integral, log_mass[n_sims], phi, alpha, delta, psi, mu[3], gene_avg;
  
  for (int g=0; g<G; g++) {
    partial_integral = 0;
    
    gene_avg = 0;
    for (int ii=0; ii<n; ii++) gene_avg += count(g,ii);
    gene_avg /= n;
    
    for (int s=0; s<n_sims; s++) {
      phi   = rnorm(1, phi_location,   phi_scale)[0]; 
      //phi   = rnorm(1, log(gene_avg),   .1)[0]; 
      alpha = rnorm(1, alpha_location, alpha_scale)[0]; // change to laplace
      delta = rnorm(1, delta_location, delta_scale)[0]; // change to laplace
      psi   = rnorm(1, psi_location,   psi_scale)[0];
      
      mu[0] = phi+alpha; // make sure these agree with definition
      mu[1] = phi-alpha;
      mu[2] = phi+delta;
      
      log_mass[s] = 0;
      for (int i=0; i<n; i++) {
        log_mass[s] += dnbinom_mu(count(g,i), 1/exp(psi), exp(mu[variety[i]]), 1);
        //Rprintf("%i %i %i\n %f %f %f %f\n %f %f %f\n %i %f\n", 
        //        g, s, i, phi, alpha, delta, psi, mu[0], mu[1], mu[2], count(g,i),dnbinom_mu(count(g,i), 1/exp(psi), exp(mu[variety[i]]), 1));
        //Rprintf("  %i %i %i %f\n", g, s, i, log_mass);
      }
    }
    integral += logsumexp(log_mass, n_sims);
  }     
         
  return integral;
}
