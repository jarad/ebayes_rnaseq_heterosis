#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

void rpadp(double *padp, NumericVector mean, NumericVector scale) {
  for (int i=0; i<4; i++) 
    padp[i] = rnorm(1, mean[i], scale[i])[0];
}

double dpadp(double *padp, NumericVector mean, NumericVector scale) {
  double sum=0;
  
  for (int i=0; i<4; i++)
    sum += R::dnorm(padp[i], mean[i], scale[i], 1);

  return sum;
}

void calculate_mu(NumericVector mu, double padp[], IntegerVector variety) {
  double tmp_mu[3];
  tmp_mu[0] = padp[0]+padp[1]; // make sure these agree with definition
  tmp_mu[1] = padp[0]-padp[1];
  tmp_mu[2] = padp[0]+padp[2];
  
  for (int i=0; i<mu.length(); i++)
    mu[i] = tmp_mu[variety[i]];
}

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
NumericVector dnbinom_mu(NumericMatrix::Column x, double size, NumericVector mu, int lg) {
  // Vectorized for mu (mean)
  NumericVector result(mu.length());
  for (int i=0; i<mu.length(); i++) 
    result[i] = dnbinom_mu(x[i], size, mu[i], lg);
  return result;
}



// [[Rcpp::export]]
double monte_carlo_integral(IntegerMatrix count, 
                            IntegerVector variety, 
                            NumericVector mean,
                            NumericVector scale,
                            int n_sims) {
  // count is the observed count (sample x gene, column major for quick access)
  // variety is the variety identifier for the columns
  // mean and scale contain the hyperparameters for the phi-alpha-delta-psi distributions
  // n_sims number of Monte Carlo simulations to use
  int G = count.ncol(); // number of genes
  int n = count.nrow(); // number of samples
  double log_n_sims = log( (double) n_sims);
  
  // argument error handling
  if (variety.size() != n) throw std::range_error("Length of variety and number of columns of count do not agree.");
//  if (variety.max()   > 3) throw std::range_error("Values in variety must not be greater than 3.")
  if (mean.size()   != 4) throw std::range_error("mean  must have length 4.");  
  if (scale.size()  != 4) throw std::range_error("scale must have length 4.");
         
  double integral = 0, // value of the integral (updated for each gene)
    log_mass[n_sims],  // log_mass for each simulated set of parameters (reused for each gene)
    padp[4];           // array containing values for phi, alpha, delta, and psi           
  
  NumericVector mu(n); // array to contain values for mu derived from padp
  
  for (int g=0; g<G; g++) {
    //Rprintf("%i: ", g); 
    
    for (int s=0; s<n_sims; s++) {
      rpadp(padp, mean, scale);        // random draws for phi, alpha, delta, and psi
      calculate_mu(mu, padp, variety); // calculate mu from phi, alpha, delta
      
      log_mass[s] = dpadp(padp, mean, scale) - // target 
                    dpadp(padp, mean, scale) + // proposal
                    dnbinom_mu(count(_,g), 1/exp(padp[3]), exp(mu), 1);  
     
//     for (int i=0; i<n; i++) {
//        log_mass[s] += dnbinom_mu(count(g,i), 1/exp(padp[3]), exp(mu[variety[i]]), 1);
//      } 
      //Rprintf("%e ", log_mass[s]);
    }
    
    integral += logsumexp(log_mass, n_sims) - log_n_sims;
    //Rprintf(", total=%f", logsumexp(log_mass, n_sims) - log_n_sims);
    //Rprintf("\n");
  }     
         
  return integral;
}
