library(edgeR)
library(plyr)
library(dplyr)
library(reshape2)

source("Laplace.R")
source("sim_heterosis_data.R")

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

tmp = 
rdply(100, {
  fit = sim_heterosis_data(G, nv=4, parameters = p)$data[,c("gene","sample","count")]  %>% 
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
}, .progress='text', .id=NULL)

saveRDS(tmp,file="savedtmp.rds")


# Calculate bias and MSE
p$type = 'truth'
tmp$type = 'estimate'
stats = rbind(p,select(tmp)) %>%
  melt(id.vars = c('gene','type'), variable.name='parameter') %>%
  group_by(gene,parameter) %>%
  summarize(truth = unique(value[type=='truth']),
            bias  = mean(value[type=='estimate'] - value[type=='truth']), 
            mse   = mean((value[type=='estimate'] - value[type=='truth'])^2))

library(ggplot2)
ggplot(stats %>% melt(id.vars=c('gene','parameter','truth'), variable.name='measure'), 
       aes(truth,value)) +
  geom_point() +
  facet_grid(measure~parameter, scales='free')

# Calculate method of moments estimators
tmp$sim = rep(1:G,each=100)
MoM = select(tmp,-type) %>%
  melt(id.vars=c('gene','sim'), variable.name='parameter') %>%
  group_by(sim,parameter) %>%
  summarize(est_mean = mean(value,na.rm=TRUE),
            est_scl   = sd(value)) %>%
  mutate(est_scl = ifelse(parameter %in% c('alpha','delta'), est_scl / sqrt(2), est_scl))

hyp.truth = data.frame(parameter = c("phi","alpha","delta","psi"),
                       tr.mean   = c(4.6,0,0,-2),
                       tr.scale  = c(1.8,.1,.01,.1))
ggplot(MoM, aes(x=est_mean))+geom_histogram()+facet_wrap(~parameter,scales="free")+
  geom_vline(data=hyp.truth, color = "red", aes(xintercept=tr.mean))

ggplot(MoM, aes(x=est_scl))+geom_histogram()+facet_wrap(~parameter,scales="free")+
  geom_vline(data=hyp.truth, color = "red", aes(xintercept=tr.scale))





# Looking at sampling distributions of estimates for single genesopar = par(mfrow=c(2,2))
i = sample(G,1)
hist(tmp$phi[tmp$gene==i]); abline(v=p$phi[p$gene==i], col='red')
hist(tmp$alpha[tmp$gene==i]); abline(v=p$alpha[p$gene==i], col='red')
hist(tmp$delta[tmp$gene==i]); abline(v=p$delta[p$gene==i], col='red')
hist(tmp$psi[tmp$gene==i]); abline(v=p$psi[p$gene==i], col='red')
par(opar)
