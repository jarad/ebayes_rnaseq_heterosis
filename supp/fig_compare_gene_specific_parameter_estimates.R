library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)
library(ggplot2)

##################################################
# Get hyperparameters
##################################################
d = read.csv("data.csv", row.names=1)
variety = factor(gsub("_[0-9]{1,2}", "", names(d)), levels=c("B73","Mo17","B73xMo17"))

# phi, alpha, delta parameterization
design = cbind(1,
               c(1,-1,0)[as.numeric(variety)],
               variety == 'B73xMo17')

# GLM fit using edgeR
fit = d %>% 
  DGEList() %>%
  calcNormFactors %>%
  estimateCommonDisp %>%
  estimateGLMTagwiseDisp(design) %>%
  glmFit(design)

# Calculate gene-specific estimates for phi, alpha, delta, and psi
hat = data.frame(phi   = fit$coefficients[,1] + mean(fit$offset[1,]),
                 alpha = -fit$coefficients[,2],
                 delta = fit$coefficients[,3],
                 psi   = log(fit$dispersion),
                 method="edgeR")
hat$gene = rownames(hat)


# Read eBayes parameter estimates
tmp = readRDS("script.rds")[,c("phi","alpha","delta","psi")]
tmp$gene = rownames(tmp)
tmp$method = 'eBayes'

# Combine estimates
combined = rbind(melt(hat, id.vars=c("method","gene")),
                 melt(tmp, id.vars=c("method","gene")))
casted = dcast(combined, gene+variable~method)

casted$variable = revalue(casted$variable, c("phi"   = "parental average",
                                             "alpha" = "half-parental difference",
                                             "delta" = "hybrid effect",
                                             "psi"   = "overdispersion"))

# Plot estimates
ggplot(casted, aes(edgeR, eBayes)) + 
  stat_binhex() + 
  geom_abline(intercept=0,slope=1) +
  facet_wrap(~variable, scales='free')
