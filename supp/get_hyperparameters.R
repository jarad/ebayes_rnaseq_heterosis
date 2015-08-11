get_hyperparameters = function(d, variety, method) {
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
  
  # Fix sd for Laplace priors
  if (method=='laplace') {
    ss$sd[ss$parameter=='alpha'] = ss$sd[ss$parameter=='alpha']/sqrt(2)
    ss$sd[ss$parameter=='delta'] = ss$sd[ss$parameter=='delta']/sqrt(2)
  }
  
  list(eta_phi   = ss$mean[ss$parameter=='phi'],
       eta_alpha = ss$mean[ss$parameter=='alpha'],
       eta_delta = ss$mean[ss$parameter=='delta'],
       eta_psi   = ss$mean[ss$parameter=='psi'],
       sigma_phi   = ss$sd[ss$parameter=='phi'],
       sigma_alpha = ss$sd[ss$parameter=='alpha'], 
       sigma_delta = ss$sd[ss$parameter=='delta'], 
       sigma_psi   = ss$sd[ss$parameter=='psi'],
       c = fit$offset[1,] - mean(fit$offset[1,]))
}
