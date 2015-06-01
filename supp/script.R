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
# This will take a long time, so parallelize if possible
if (parallel <- require(doMC)) {
  registerDoMC()
} else if (parallel <- require(doParallel)) {
  cl = makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
}

single_gene_analysis = function(.data) {
  diverge = TRUE
  attempt = 1
  
  
  while (!divere) {
    r = sampling(model, 
                 c(list(S       = nrow(.data), 
                        count   = .data$value,
                        variety = as.numeric(factor(gsub("_[1-4]", "", .data$sampleID), levels=c("B73","Mo17","B73xMo17")))),
                   hyperparameters),
                 iter = 2000*attempt,
                 thin = attempt)
    # Check PSRF for (lack of) convergence
    diverge = any()
    attempt = attempt + 1
  }
  
  
  
  e = extract(r)
  prob_HPH = mean(abs(e$delta) >  abs(e$alpha))
  prob_LPH = mean(abs(e$delta) < -abs(e$alpha))
}


analysis = d %>% 
  mutate(gene = 1:nrow(d)) %>%
  melt(id.vars = 'gene', variable.name = 'sampleID') %>%
  subset(gene<11) %>%
  group_by(gene) 
  
