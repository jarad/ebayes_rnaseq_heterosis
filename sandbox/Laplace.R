#' Laplace (double exponential distribution)
#' 
#' Density, distribution function, quantile function, and random generation for the Laplace 
#' (double exponential) distribution with mean equal to mean and scale equal to sd.
#' 
#' @usage dlaplace(x, mean = 0, scale = 1, log=FALSE)
#' @usage rlaplace(n, mean = 0, scale = 1)
#' 
#' @param x vector of quantiles.
#' @param n number of observations. If length(x) > 1, the length is taken to be the number required.
#' @param mean vector of means.
#' @param scale vector of scales.
#' 
#' @details The Laplace distribution has density ...
#' 
#' @return dlaplace gives the density, plaplace gives the distribution function, 
#'   qlaplace gives the quantile function, and rlaplace generates random deviates. 
#'   
#'   The length of the result is determined by n for rlaplace, and is the maximum 
#'   of the lengths of the numerical arguments for the other functions.
#'   
#'   The numerical arguments other than n are recycled to the length of the result. 
#'   Only the first elements of the logical arguments are used.
#'   
#' @examples dlaplace(0)
#' 
dlaplace = function(x, mean=0, scale=1, log=FALSE) {
  logd = -abs(x-mean)/scale - log(2*scale)
  if (log) return(logd)
  return(exp(logd))
}

rlaplace = function(n, mean=0, scale=1) {
  u = runif(n)-.5
  mean - scale * sign(u)*log(1-2*abs(u))
}
