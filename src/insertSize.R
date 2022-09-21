#!Rscript

## this is for plotting insertSize

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 4) stop("Wrong number of input parameters: 'path sampleName insFile outPDF'")

path <- inputpar[1]
sampleName <- inputpar[2]
insFile <- inputpar[3]
outPDF <- inputpar[4]

setwd(path)

insert.size = read.table(gzfile(insFile))
dd = density(insert.size[,1])
med = median(insert.size[,1])
pdf(file = outPDF, width=6, height=6)
plot(dd, main=paste(sampleName, "Insert Size", sep=" "), xlim=c(0,500))
points(med, max(dd$y), pch=20, col=rgb(1,0,0,1/2), cex=1.6)
text(med, quantile(dd$y, prob=seq(0,1,0.01))["99%"], labels=paste("median=", round(med,1), sep=""))
dev.off()

bb = boxplot(insert.size[,1], plot=F)
save(bb,dd,file=paste(path, sampleName, ".ins.rda", sep=""))
