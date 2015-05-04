d = read.table("../shiny/rna_seq.txt", header=TRUE)

y = d[,grep('total',names(d))]

filter = apply(y, 1, function(x) length(x[x>5])>=2)

filtered = y[filter,]
rownames(filtered) = d$GeneID[filter]
spikes = rownames(filtered)[grep("AC", rownames(filtered))]



x <- as.factor(rep(c("B","M","BM","MB"),4))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


set1 <- RUVg(set, spikes, k=1)
pData(set1)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)
