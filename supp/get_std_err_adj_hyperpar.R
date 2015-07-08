get_std_err_adj_hyperpar = function(d, variety) {
  
  require(MASS)
  
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
  
  theta = 1/fit$dispersion
  offset = fit$offset[1,] - mean(fit$offset[1,])
  
  med_se = ldply(1:nrow(d), function(gene){
    wts = glm.fit(x=design,
                 y=as.numeric(d[gene,]),
                 offset = offset,
                 family=negative.binomial(theta=theta[gene]))$weights
    
    se = sqrt(diag(solve(t(design) %*% diag(wts) %*% design))) #std errors
    se = as.data.frame(t(se))
    names(se) = c('phi','alpha','delta')
    se
    }) %>%
    summarise(phi = median(phi),
              alpha = median(alpha),
              delta = median(delta))
  
  ss = hat %>%
    melt(id.vars='gene', variable.name = 'parameter') %>%
    group_by(parameter) %>%
    summarize(mean=mean(value), sd=sd(value))
  
  list(eta_phi   = ss$mean[ss$parameter=='phi'],
       eta_alpha = ss$mean[ss$parameter=='alpha'],
       eta_delta = ss$mean[ss$parameter=='delta'],
       eta_psi   = ss$mean[ss$parameter=='psi'],
       sigma_phi   = sqrt(ss$sd[ss$parameter=='phi']^2 - med_se$phi^2),
       sigma_alpha = sqrt(ss$sd[ss$parameter=='alpha']^2 - med_se$alpha^2)/sqrt(2), # Fix sd for Laplace priors
       sigma_delta = sqrt(ss$sd[ss$parameter=='delta']^2 - med_se$delta^2)/sqrt(2), # Fix sd for Laplace priors
       sigma_psi   = ss$sd[ss$parameter=='psi'],
       c = offset)
}
