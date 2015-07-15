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
