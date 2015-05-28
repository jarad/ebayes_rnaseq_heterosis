library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(MASS) #?
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

# Fix standard deviations for Laplace priors
ss$sd[2:3] = ss$sd[2:3]/sqrt(2)


######################################
# Run individual gene analyses
######################################
# This will take a long time, so parallelize if possible
if (parallel <- require(doMC)) {
  registerDoMC()
} 

analysis = d %>% 
  mutate(gene = 1:nrow(d)) %>%
  melt(id.vars = 'gene', variable.name = 'sampleID') %>%
  group_by(gene) 
  
