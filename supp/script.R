library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)

# Compile stan model
model = stan_model("model.stan")

d = read.csv("data.csv", row.names=1)
variety = factor(gsub("_[1-4]", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# trim genes for average count greater than 1
d = d[which(rowMeans(d)>1),] 

#######################################
# Find hyperparameter estimates
#######################################
# phi, alpha, delta parameterization
design = cbind(1,
               c(1,-1,0)[as.numeric(variety)],
               variety == 'B73xMo17')

# GLM fit using edgeR
fit = d %>% 
  DGEList() %>%
  calcNormFactors %>%
  estimateCommonDisp(design) %>%
  estimateGLMTagwiseDisp(design) %>%
  glmFit(design)

# Calculate gene-specific estimates for phi, alpha, and delta
hat = data.frame(gene = 1:length(fit$dispersion),
                 phi   = fit$coefficients[,1] + mean(fit$offset[1,]),
                 alpha = fit$coefficients[,2],
                 delta = fit$coefficients[,3],
                 psi   = log(fit$dispersion))

ss = hat %>%
  melt(id.vars='gene', variable.name = 'parameter') %>%
  group_by(parameter) %>%
  summarize(mean=mean(value), sd=sd(value))

hyperparameters = list(eta_phi   = ss$mean[ss$parameter=='phi'],
                       eta_alpha = ss$mean[ss$parameter=='alpha'],
                       eta_delta = ss$mean[ss$parameter=='delta'],
                       eta_psi   = ss$mean[ss$parameter=='psi'],
                       sigma_phi   = ss$sd[ss$parameter=='phi'],
                       sigma_alpha = ss$sd[ss$parameter=='alpha']/sqrt(2), # Fix sd for Laplace priors
                       sigma_delta = ss$sd[ss$parameter=='delta']/sqrt(2), # Fix sd for Laplace priors
                       sigma_psi   = ss$sd[ss$parameter=='psi'],
                       c = fit$offset[1,] - mean(fit$offset[1,]))


######################################
# Run individual gene analyses
######################################


single_gene_analysis = function(counts) {
  diverge = TRUE
  attempt = 1
  
  
  while (diverge) {
    r = sampling(model, 
                 data = c(list(S       = length(counts), 
                               count   = as.numeric(counts),
                               variety = as.numeric(factor(gsub("_[1-4]", "", colnames(counts)), levels=c("B73","Mo17","B73xMo17")))),
                          hyperparameters),
                 pars = c("phi","alpha","delta","psi","LPH","HPH","effectiveSize2"),
                 iter = 2000*attempt,
                 thin = attempt)
    
    # Check PSRF for (lack of) convergence
    s = summary(r)$summary
    diverge = any(s[,"n_eff"] < 1000)
    attempt = attempt + 1
  }
  
  data.frame(
    prob_LPH = s[rownames(s) == "LPH","mean"],
    prob_HPH = s[rownames(s) == "HPH","mean"],
    effective = s[rownames(s) == "effectiveSize2","mean"])
}



# This will take a long time, so parallelize if possible
if (parallel <- require(doMC)) {
  registerDoMC()
} else if (parallel <- require(doParallel)) {
  cl = makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
}

wh = 1:2 # which genes to run the analysis on
analysis = adply(d[wh,],
                 1,
                 function(x) single_gene_analysis(x),
                 .id = 'gene',
                 .parallel = parallel,
                 .paropts = list(.export=c('single_gene_analysis','model','hyperparameters'), .packages='rstan'))

