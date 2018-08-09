#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
print(args)
myfile = 'stdin.txt'

 d <- read.table(myfile, sep = "\t", header = F)
if (nrow(d) > 0){
    pvalues <- vector(mode="double", length=dim(d)[1])
    oddratio <- vector(mode="double", length=dim(d)[1])
    for( i in 1:dim(d)[1] ) {
        h <- fisher.test(matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2))
        pvalues[i] <- round(h$p.value, 5)
        oddratio[i] <- round(h$estimate, 5)
    }
    write.table(data.frame(d[,1:20], pvalues, oddratio, d[,21:dim(d)[2]]), file="", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
}
