d = read.table("rna_seq.txt", header=TRUE)

d = d[,c(1, grep("total", names(d)))]

names(d)[-1] = paste(c("B73","Mo17","B73xMo17","Mo17xB73"), 
                     rep(1:4, each=4), 
                     sep="_")

# Order columns
B_cols = c(2,6,10,14)
d = d[,c(1, B_cols, B_cols+1, B_cols+2, B_cols+3)]

# Assuming we are using the BM cross
d = d[,1:13]

# Write file
write.csv(d, file="../supp/data.csv", row.names=FALSE)
