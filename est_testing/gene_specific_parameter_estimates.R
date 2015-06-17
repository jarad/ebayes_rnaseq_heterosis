library(edgeR)
library(plyr)
library(dplyr) #must load dplyr second
library(reshape2)

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

set.seed(101669)

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
  hat
}, .progress='text', .id=NULL)

saveRDS(tmp,file="savedtmp.rds")
# tmp = readRDS("savedtmp.rds")

#Include standard errors for parameters for a single simulation
set.seed(691152)
tmp2 = 
  rdply(20, {
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
    hat$se_phi   = abs(hat$phi)/sqrt(glmLRT(fit,coef=1)$table[,3]+.0001)
    hat$se_alpha = abs(hat$alpha)/(sqrt(glmLRT(fit,coef=2)$table[,3])+.0001)
    hat$se_delta = abs(hat$delta)/(sqrt(glmLRT(fit,coef=3)$table[,3])+.0001)
    hat
}, .progress='text', .id=NULL)

tmp2$sim = rep(1:20,each=G)
adjscale = select(tmp2,-gene) %>%
  ddply(.(sim),summarise,
    phi   = sqrt(max(c(var(phi)   - median(se_phi)^2, 0))),
    alpha = sqrt(max(c(var(alpha) - median(se_phi,na.rm=T)^2, 0)))/sqrt(2),
    delta = sqrt(max(c(var(delta) - median(se_delta,na.rm=T)^2, 0)))/sqrt(2))

hyp.truth = data.frame(parameter = c("phi","alpha","delta"),
                       tr.mean   = c(4.6,0,0),
                       tr.scale  = c(1.8,.1,.01))

select(adjscale,-sim) %>% melt(variable.name="parameter") %>%
  ggplot(aes(x=value))+geom_histogram()+facet_wrap(~parameter,scales="free")+
  geom_vline(data=hyp.truth,color="red",aes(xintercept=tr.scale))


# Calculate bias and MSE
p$type = 'truth'
tmp$type = 'estimate'
stats = rbind(p,tmp) %>%
  melt(id.vars = c('gene','type'), variable.name='parameter') %>%
  group_by(gene,parameter) %>%
  summarize(truth = unique(value[type=='truth']),
            bias  = mean(value[type=='estimate'] - value[type=='truth']), 
            mse   = mean((value[type=='estimate'] - value[type=='truth'])^2))


library(ggplot2)
ggplot(stats %>% melt(id.vars=c('gene','parameter','truth'), variable.name='measure'), 
       aes(truth,value)) +
  geom_hex(bins=50) +
  facet_grid(measure~parameter, scales='free') +
  stat_smooth(method="loess",color="red")
ggsave("bias_mse.pdf")
# Calculate method of moments estimators
tmp$sim = rep(1:100,each=G)
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

# Looking at parameter bias, mse conditional on true phi
  filter(stats, parameter != "phi") %>%
  select(-mse) %>%              
  dcast(gene ~ parameter) %>%
  cbind(phi=p$phi) %>%
  melt(id.vars=c("gene","phi"),value.name="bias") %>%
  ggplot(aes(x=phi,y=bias)) +
    geom_hex(bins=50)+ 
    facet_wrap(~variable,scales="free")+
    stat_smooth(method="loess",color="red")
  
  filter(stats, parameter != "phi") %>%
    select(-bias) %>%              
    dcast(gene ~ parameter) %>%
    cbind(phi=p$phi) %>%
    melt(id.vars=c("gene","phi"),value.name="mse") %>%
    ggplot(aes(x=phi,y=mse)) +
    geom_hex(bins=50)+ 
    facet_wrap(~variable,scales="free")+
    stat_smooth(method="loess",color="red")

# Looking at sampling distributions of estimates for single genesopar = par(mfrow=c(2,2))
i = sample(G,1)
hist(tmp$phi[tmp$gene==i]); abline(v=p$phi[p$gene==i], col='red')
hist(tmp$alpha[tmp$gene==i]); abline(v=p$alpha[p$gene==i], col='red')
hist(tmp$delta[tmp$gene==i]); abline(v=p$delta[p$gene==i], col='red')
hist(tmp$psi[tmp$gene==i]); abline(v=p$psi[p$gene==i], col='red')
par(opar)
