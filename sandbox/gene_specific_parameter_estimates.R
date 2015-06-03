library(edgeR)
library(plyr)
library(dplyr)
library(reshape2)

source("Laplace.R")
source("sim_heterosis_data.R")

# 
G = 250
d = sim_heterosis_data(G, nv=4, verbose=1)
p = d$parameters # Get gene specific parameter draws

for (i in 1:nrow(p)) p[,2:5] = p[84,2:5]

# phi, alpha, delta parameterization
variety = d$data$variety[d$data$gene==1]
design = cbind(1,
               c(1,-1,0)[variety],
               variety == 3)

tmp = 
rdply(10, {
  fit = sim_heterosis_data(G, parameters = p)$data[,c("gene","sample","count")]  %>% 
    dcast(formula = gene~sample, value.var='count') %>% # Assumes columns are in variety order
    select(-gene) %>%
    as.matrix %>%
    DGEList %>%
    calcNormFactors %>%
    estimateCommonDisp %>%
    estimateGLMTagwiseDisp(design) %>%
    glmFit(design)
  
  hat = data.frame(gene = 1:length(fit$dispersion),
             phi   = fit$coefficients[,1] + mean(fit$offset[1,]),
             alpha = fit$coefficients[,2],
             delta = fit$coefficients[,3],
             psi   = log(fit$dispersion))
}, .progress='text')


