
edgeR_est <- function(counts, group, geneid){

  require(edgeR)
  row.names(counts) <- geneid
  countsDGE         <- DGEList(counts)
  countsDGE         <- calcNormFactors(countsDGE)
  design            <- model.matrix(~ 0 + as.factor(group))
  countsDGE         <- estimateGLMCommonDisp(countsDGE, design)
  countsDGE         <- estimateGLMTagwiseDisp(countsDGE, design)
  fit               <- glmFit(countsDGE, design)
  offset            = fit$offset                                 # offset = log(lib.size * norm.factors)
  mean.off          = mean(offset[1,])
  offset.new        = offset - mean.off                          # offset is relative size in our model
  mu                = fit$coefficients + mean.off
  lin.comb          = matrix(c( .5, .5, 0,
                               -.5, .5, 0,
                               -.5,-.5, 1),3,3,byrow=T )
  B   = lin.comb %*% t(mu)
  psi = log(fit$dispersion)


  return(list(parameters = data.frame(phi = B[1,], alpha = B[2,], delta = B[3,], psi = psi),
              hyperparameters = data.frame(parameter = c("phi","alpha","delta","psi"),
                                           location  = c(rowMeans(B), mean(psi)),
                                           scale     = c(sd(B[1,]), sd(B[2,])/sqrt(2), sd(B[3,])/sqrt(2), sd(psi))),
              c = offset.new[1,],
              sigma_c = sd(offset.new[1,])
              ))
}

