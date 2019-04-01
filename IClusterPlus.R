if(!require(iClusterPlus)){
  source("https://bioconductor.org/biocLite.R")
  biocLite('iClusterPlus')
}
if(!require(GenomicRanges)){
  source("https://bioconductor.org/biocLite.R")
  biocLite('GenomicRanges')
}
#biocLite('GenomicRanges')
library(iClusterPlus)
library(GenomicRanges)
load('iClusterPlus.preData.RData')
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=12,dt1=mtrix.mut2,dt2=region.cnv,dt3=t(gene.exp),
                             dt4=t(methy),
          type=c("binomial","gaussian","gaussian"),K=k,n.lambda=2129,
          scale.lambda=c(1,1,1),maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}