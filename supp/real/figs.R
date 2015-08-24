library(ggplot2)
library(gridExtra)

d = readRDS("script.rds")

##########################################################################################
# Volcano plot
##########################################################################################
ggplot(d, aes(effectiveSize, pmax(prob_LPH, prob_HPH))) + 
  stat_binhex(color="white",bin=100) +
  theme_bw() + 
  scale_fill_gradientn(trans='log',breaks=c(1,10,100,1000), colours=c("gray","black")) +
  labs(x="Effect size", y="Maximum heterosis probabilities")

ggsave("../include/volcano.pdf", width=8, height=8)
  

##########################################################################################
# Hyperparameter plot
##########################################################################################
parms = c("phi","alpha","delta","psi")
estimates = d[,parms]
plots = list()

n = ncol(estimates)
for (i in 1:n) {
  for (j in 1:n) {
    p = (i-1)*n+j # plot number
    if (i>j) 
      plots[[p]] = grid.rect(gp=gpar(col="white"))
    if (i==j) {
      plots[[p]] = ggplot(estimates, aes_string(parms[i])) + 
        geom_histogram(fill='gray') +
        theme_bw() + 
        labs(x=parse(text=parms[i]))
    }
    if (i<j) {   
      plots[[p]] = ggplot(estimates, aes_string(parms[j], parms[i])) + 
        stat_binhex(bin=100) +
        theme_bw() + 
        scale_fill_gradientn(trans='log',breaks=c(1,10,100,1000), guide=FALSE, colours=c("gray","black")) + 
        labs(x=parse(text=parms[j]), y=parse(text=parms[i]))
    }
  }
}

do.call("grid.arrange", c(plots, ncol=4))

ggsave("../include/estimates.pdf", do.call("arrangeGrob", c(plots, ncol=4)), width=8, height=8)



##########################################################################################
# Comparison of gene specific parameter estimates
##########################################################################################

library(dplyr)
library(reshape2)
library(edgeR)

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
                 alpha = fit$coefficients[,2], 
                 delta = fit$coefficients[,3],
                 psi   = log(fit$dispersion),
                 method="edgeR")
hat$gene = rownames(hat)


# Read eBayes parameter estimates
tmp = readRDS("script.rds")[,c("phi","alpha","delta","psi")]
tmp$gene = rownames(tmp)
tmp$method = 'eBayes'
tmp$alpha = -tmp$alpha # Fixed in the stan model files, but the rds file has swapped parents

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
  facet_wrap(~variable, scales='free') +
  theme_bw() + 
  coord_fixed() +
  scale_fill_gradientn(trans='log',breaks=c(1,10,100,1000), colours=c("gray","black"))


ggsave(file="../include/gene_specific_estimates.pdf", width=8, height=8)



