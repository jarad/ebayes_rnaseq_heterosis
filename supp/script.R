library(rstan)
library(plyr)
library(dplyr)
library(reshape2)
library(MASS)

# Compile stan model
model = stan_model("model.stan")

d = read.csv("data_matrix.csv", row.names=1)

# trim genes for average count 
d = d[which(rowMeans(d)>0),] 

# use edgeR 
library(edgeR)

# phi, alpha, delta parameterization
design = cbind(1,
               rep(c(1,-1,0), each=4),
               rep(c(0, 0,1), each=4))

# GLM fit
fit = d %>% 
  DGEList() %>%
  calcNormFactors %>%
  estimateCommonDisp(design) %>%
  estimateGLMTagwiseDisp(design) %>%
  glmFit(design)

