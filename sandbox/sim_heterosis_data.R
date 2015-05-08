#' Simulate heterosis data
#' 
#' Simulates heterosis data from either hyperparameters (and thus simulating gene-specific parameters)
#' or gene-specific parameters.
#' 
#' @param G integer. The number of genes to simulate.
#' @param nv integer (or integer vector). indicating the number of samples per genotype
#' @param parameters. a list containing either hyperparameters or (gene-specific) parameters
#' @param distributions. a list containing the distribution to use for the gene-specific parameters
#' 
#' @return list containing the hyperparameters, gene-specific parameters, and data
#' 

sim_heterosis_data = function(G=10, nv=4, parameters=NULL, distributions=NULL) {
  require(plyr)
  
  if (length(nv==1)) nv = rep(nv,3)
  
  if (is.null(parameters)) 
    hyperparameters = data.frame(parameter      = c("phi","alpha","delta","psi"),
                                 location = c(4.6,0,0,-2),
                                 scale    = c(1.8,.1,.1,.1))
  
  # Simulate gene-specific parameters
  if (is.null(distributions)) 
    parameters = rdply(G, {
      with(hyperparameters, 
           data.frame(phi   = rnorm   (1, location[parameter == "phi"  ], scale[parameter == "phi"  ]),
                      alpha = rlaplace(1, location[parameter == "alpha"], scale[parameter == "alpha"]),
                      delta = rlaplace(1, location[parameter == "delta"], scale[parameter == "delta"]),
                      psi   = rnorm   (1, location[parameter == "psi"  ], scale[parameter == "psi"  ])))
    }, .id="gene")
  
  data = ddply(data.frame(gene=1:G), .(gene), function(x) {
    g = as.numeric(x$gene)
    mutate(data.frame(variety = rep(1:3, times = nv)),
           sample = 1:sum(nv),
           eta = with(parameters, c(phi[g]+alpha[g], phi[g]-alpha[g], phi[g]+delta[g]))[variety],
           count = rnbinom(length(eta), size = 1/exp(parameters$psi[g]), mu = exp(eta)))
  }, .inform=T)
  
  return(list(data            = data, 
              parameters      = parameters,
              hyperparameters = hyperparameters))
}
