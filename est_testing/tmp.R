library(MASS)
library(plyr)
library(dplyr)
library(reshape2)
library(edgeR)
library(stats)

parallel=FALSE
if(require(doMC)){
  parallel=TRUE
}


source("../sandbox/Laplace.R")
source("../sandbox/sim_heterosis_data.R")

set.seed(101469)

# 
G = 2500
d = sim_heterosis_data(G, nv=4, verbose=1)
p = d$parameters # Get gene specific parameter draws

# phi, alpha, delta parameterization
variety = d$data$variety[d$data$gene==1]
design = cbind(1,
               c(1,-1,0)[variety],
               variety == 3)


set.seed(691152)
res = rdply(100,{
  sim = sim_heterosis_data(G, nv=4, parameters = p)
  fit = dcast(sim$data[,c("gene","sample","count")],formula=gene~sample, value.var='count') %>%
    select(-gene) %>%
    as.matrix %>%
    DGEList %>%
    calcNormFactors %>%
    estimateCommonDisp %>%
    estimateGLMTagwiseDisp(design)%>%
    glmFit(design)
  
  hat = data.frame(gene = 1:length(fit$dispersion),
                   phi   = fit$coefficients[,1] + mean(fit$offset[1,]),
                   alpha = fit$coefficients[,2],
                   delta = fit$coefficients[,3],
                   psi   = log(fit$dispersion))
  
  theta = 1/fit$dispersion
  SE = ldply(1:G, function(i){
    one_g = filter(sim$data,gene==i)
    X=matrix(c(rep(c(1,1,1),each=4),rep(c(-1,1,0),each=4),rep(c(0,0,1),each=4)),12,3)
    ffit = glm.fit(x=X,
                  y=one_g$count,offset = fit$offset[1,]-mean(fit$offset[1,]),
                  family=negative.binomial(theta=theta[i]))
    se = sqrt(diag(solve(t(X)%*%diag(ffit$weights)%*%X))) #std errors
  })
  cbind(hat,SE)
}, .progress='text', .parallel=parallel, .id=NULL)



res$sim = rep(1:50,each=G)
res$gene = rep(1:G,times=50)
emp_sd = ddply(res, .(gene), summarise,
      sd_ahat = sd(X2),
      sd_dhat = sd(X3))

g = sample(G, 1)
filter(res, gene == g) %>%
  ggplot(aes(x = res$V3)) + geom_histogram(binwidth=.02) +
  geom_vline(xintercept=emp_sd$sd_dhat[g]) +
  scale_x_log10() + xlim(c(0,1))

emp_tr_se = melt(res,id.vars="sim",value.name = "mean") %>%
  ddply(.(sim,par),summarise,
        sd())

round(cbind(sd(fit$coefficients[,3]),quantile(SE$V3,0:10*.1)),3)



