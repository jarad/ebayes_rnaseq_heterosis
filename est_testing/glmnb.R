library(MASS)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)
library(stats)

parallel=FALSE
if(require(doMC)){
  registerDoMC(8)
  parallel=TRUE
}


source("../sandbox/Laplace.R")
source("../sandbox/sim_heterosis_data.R")

set.seed(101469)

# 
G = 2500
hyp.truth = data.frame(parameter = c("phi","alpha","delta","psi"),
                       tr.mean   = c(4.6,0,0,-2),
                       tr.scale  = c(1.8,.1,.01,.1))
d = sim_heterosis_data(G, nv=4, hyperparameters = hyp.truth, verbose=1)
p = d$parameters # Get gene specific parameter draws

# phi, alpha, delta parameterization
variety = d$data$variety[d$data$gene==1]
design = cbind(1,
               c(1,-1,0)[variety],
               variety == 3)

nb.robust = function(data, model){
  tryCatch(glm.nb(count ~ model, data=data),
              warning = function(w){
                glm(count ~ model, family = poisson(link="log"), data=data)
              })
}

extract.fit = function(f){
  coef = c(f$coef)
  se   = c(sqrt(diag(vcov(f))))
  method = as.character(f$call[[1]])
  
  if(method =="glm") {
    coef = c(coef,Inf,0)
    se   = c(se, Inf, 0)
  } else {
    coef = c(coef,f$theta,-log(f$theta))
    se   = c(se,f$SE.theta, 1/f$theta * f$SE.theta) #delta method
  }
  
  data.frame(par = c("phi","alpha","delta","theta","psi"), coef, se, method)
}

valid.res = function(fit){
  tryCatch(extract.fit(fit),
    error = function(e){
      data.frame(par=NA,coef=NA,se=NA,method="error")
    }
  )
}

set.seed(691152)
res = ldply(1:100, function(r){

  sim = sim_heterosis_data(G, nv=4, parameters = p)

  fit = ddply(sim$data, .(gene), function(g){
          f = tryCatch(nb.robust(data = g, model = design[,2:3]),
                       error = function(e){NULL})
          return(valid.res(f))
    })
  fit$sim = r
  return(fit)
}, .progress = 'text', .parallel = parallel)



