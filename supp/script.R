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
dge = d %>% 
  DGEList() %>%
  calcNormFactors %>%
  estimateCommonDisp(design) %>%
  estimateGLMTagwiseDisp(design) %>%
  glmFit(design)




# long format
m = melt(d, 
         id.var = "GeneID",
         variable.name="sample", 
         value.name="count")
m$variety = factor(gsub("_[1-4]", 
                        "", 
                        as.character(m$sample)),
                   levels=c("B73","Mo17","B73xMo17"))

ind = ddply(m, .(GeneID), function(x) {
  mo = glm.nb(count+1~variety, data=x)
  as.data.frame(coef(mo))
}, .inform=TRUE)

r = glm(count+1~variety,data=test, family="poisson")

tmp = 
m %>% 
  group_by(GeneID, variety) %>%
  summarize(mean=mean(count)) %>%
  dcast(GeneID~variety, value.var="mean")

zeros = 
tmp %>% 
  rowwise() %>%
  summarize(any = B73==0 | Mo17 == 0 | B73xMo17 == 0)
mean(zeros$any)
sum(zeros$any)


tmp[zeros$any,] 
  


rm = rowMeans(d[,-1])

mean(rowM)