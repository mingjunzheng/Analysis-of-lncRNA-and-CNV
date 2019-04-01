#source("https://bioconductor.org/biocLite.R")
#biocLite('GenomicRanges')
setwd('D:/work/MG20181215_lnc_cnv_ov')

library(iClusterPlus)
library(GenomicRanges)
library(lattice)
gene.exp.file='pcg.exp.txt'
lnc.exp.file='lnc.exp.txt'
mut.file='TCGA.OV.mutect.b22b85eb-2ca8-4c9f-a1cd-b77caab999bd.DR-10.0.somatic.maf'
cnv.file='Copy Number Variation.nocnv.merge.txt'
methy.file='HumanMethylation27'
methy.inCpG.file='CpG.Anno.clear.txt'
methy.outCpG.file='CpG_CrossMap.txt'
rnaseq.count.file='HTSeq - Counts.merge.txt'
gene.type.file='GeneTag.txt'
clinical.file='Merge_mongodb.txt'
###############Clinical#################
clinical=read.csv(clinical.file,sep = '\t',check.names = F,stringsAsFactors = F)
clinical=clinical[clinical$cancer_code=='OV',]
dim(clinical)

clinical$sample_name
clinical$os
clinical$os_status
library('survival')
library('survivalROC')
coxFun=function(dat){
  colnames(dat)=c('time','status','gene')
  fmla <- as.formula("Surv(time, status) ~gene")
  cox <- coxph(fmla, data =dat)
  p <- summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
plotKMCox=function(dat){
  colnames(dat)=c('time','status','groups')
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  print((sdf))
  p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  colKm=rainbow(length(sf$strata))
  plot(sf, mark.time = TRUE,col=colKm,xlab=paste("Survival time in day","\np=",round(p,5)),ylab = "Survival probabilities",main="Method Kaplan Meier")
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(N=',sdf$n,')'), col = colKm,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE,cex = 0.8)
  return(p)
}

sType='-01'
###############Pre gene#################
gene.exp=read.csv(file = gene.exp.file,sep = '\t',row.names = 1,check.names = F,stringsAsFactors = F)
dim(gene.exp)
comm.smps=intersect(paste0(clinical$sample_name,sType),
  colnames(gene.exp))

gene.cox=t(apply(gene.exp, 1,function(x){
  vl=as.numeric(x[match(comm.smps,colnames(gene.exp))])
  tm=clinical$os[match(comm.smps,paste0(clinical$sample_name,sType))]
  ev=clinical$os_status[match(comm.smps,paste0(clinical$sample_name,sType))]
  dat=data.frame(tm,ev,vl)[which(tm>30&!is.na(ev)&!is.na(vl)),]
  #print(dim(dat))
  if(nrow(dat)>0&sd(dat[,3])>0){
    return(coxFun(dat))
  }else{
    return(c(1,1,1,1))
  }
}))
head(gene.cox)
gene.exp=gene.exp[which(gene.cox[,1]<0.05),]
  

#length(which(apply(gene.exp, 1, mad)>quantile(apply(gene.exp, 1, mad))['75%']&
#        apply(gene.exp, 1, sd)>quantile(apply(gene.exp, 1, sd))['75%']))
#gene.exp=gene.exp[which(apply(gene.exp, 1, mad)>quantile(apply(gene.exp, 1, mad))['75%']&
#                          apply(gene.exp, 1, sd)>quantile(apply(gene.exp, 1, sd))['75%']),]


###############Pre CNV###############
cnv=read.csv(file = cnv.file,sep = '\t',check.names = F,stringsAsFactors = F)
write.table(cnv[grep('-01$',cnv[,1]),],file = 'gistic.cnv.txt',sep = '\t',quote = F,row.names = F)
region.cnv=CNregions(seg=cnv,epsilon=0,adaptive=FALSE,
                 frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
length(unique(cnv[,1]))
region.cnv[,1]

comm.smps=intersect(paste0(clinical$sample_name,sType),
                    row.names(region.cnv))
region.cnv.cox=t(apply(region.cnv, 2,function(x){
  vl=as.numeric(x[match(comm.smps,row.names(region.cnv))])
  tm=clinical$os[match(comm.smps,paste0(clinical$sample_name,sType))]
  ev=clinical$os_status[match(comm.smps,paste0(clinical$sample_name,sType))]
  dat=data.frame(tm,ev,vl)[which(tm>30&!is.na(ev)&!is.na(vl)),]
  #print(dim(dat))
  if(nrow(dat)>0&sd(dat[,3])>0){
    return(coxFun(dat))
  }else{
    return(c(1,1,1,1))
  }
}))
dim(region.cnv.cox)
region.cnv=region.cnv[,which(region.cnv.cox[,1]<0.05)]

###############Pre Mut######################

#####################Methylation################################
methy=read.csv(file = methy.file,sep = '\t',check.names = F,stringsAsFactors = F,row.names = 1)
dim(methy)
methy=methy[which(apply(methy, 1, function(x){return(sum(is.na(x)))})==0),]
median(apply(methy, 1, function(x){return(sd(x))}))
###450k###
methy.in=read.csv(file = methy.inCpG.file,sep = '\t',check.names = F,stringsAsFactors = F)
methy.out=read.csv(file = methy.outCpG.file,sep = '\t',check.names = F,stringsAsFactors = F)
lst.cpgs=unique(setdiff(methy.in[,1],methy.out[,1]))
methy=methy[match(intersect(lst.cpgs,row.names(methy)),row.names(methy)),]
#methy=methy[which(apply(methy, 1, function(x){return(sd(x))})>median(apply(methy, 1, function(x){return(sd(x))}))),]##Cut CpGs
dim(methy)
#methy=methy[which(apply(methy, 1, mad)>quantile(apply(methy, 1, mad))['75%']&
#                          apply(methy, 1, sd)>quantile(apply(methy, 1, sd))['75%']),]
#head(methy)
#dim(methy)
comm.smps=intersect(paste0(clinical$sample_name,sType),
                    colnames(methy))
methy.cox=t(apply(methy, 1,function(x){
  vl=as.numeric(x[match(comm.smps,colnames(methy))])
  tm=clinical$os[match(comm.smps,paste0(clinical$sample_name,sType))]
  ev=clinical$os_status[match(comm.smps,paste0(clinical$sample_name,sType))]
  dat=data.frame(tm,ev,vl)[which(tm>30&!is.na(ev)&!is.na(vl)),]
  #print(dim(dat))
  if(nrow(dat)>0&sd(dat[,3])>0){
    return(coxFun(dat))
  }else{
    return(c(1,1,1,1))
  }
}))
save(methy.cox,file='step1.methy.cox.RData')
methy=methy[which(methy.cox[,1]<0.05),]
#sum(methy.cox[,1]<0.05)
####################Match Samples###############################
comman.samples=colnames(gene.exp)##
comman.samples=intersect(comman.samples,row.names(region.cnv))
#comman.samples=intersect(comman.samples,row.names(mtrix.mut2))
comman.samples=intersect(comman.samples,colnames(methy))
comman.samples=comman.samples[grep('-01$',comman.samples)]
print(paste0('Comman Sample:',length(comman.samples)))
gene.exp=gene.exp[,match(comman.samples,colnames(gene.exp))]
#save(gene.exp,file='gene.exp.RData')
region.cnv=region.cnv[match(comman.samples,row.names(region.cnv)),]
#save(region.cnv,file='region.cnv.RData')
#mtrix.mut2=mtrix.mut2[match(comman.samples,row.names(mtrix.mut2)),]
#save(mtrix.mut2,file='mtrix.mut2.RData')
methy=methy[,match(comman.samples,colnames(methy))]
#save(methy,file='methy.RData')
dim(gene.exp)
dim(region.cnv)
dim(methy)

save(gene.exp,region.cnv,methy,file = 'iClusterPlus.preData.RData')
###############################################################
#load('gene.exp.RData')
#load('region.cnv.RData')
#load('mtrix.mut2.RData')
#load('methy.RData')

#intersect(colnames(mtrix.mut2),row.names(lnc.exp))
#table(mtrix.mut[,which(colnames(mtrix.mut)=='MALAT1')])
#table(mtrix.mut[,which(colnames(mtrix.mut)=='XIST')])

#source("https://bioconductor.org/biocLite.R")
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
dim(region.cnv)
dim(gene.exp)
dim(methy)

fit.single=iClusterPlus(dt1=region.cnv,dt2=t(gene.exp),dt3=t(methy),
                        type=c("gaussian","gaussian","gaussian"),
                        lambda=c(0.04,0.61,0.90),K=5,maxiter=10)
save(fit.single,file='fit.single.RData')

load('fit.single.RData')
write.table(cbind(Cluster=names(table(fit.single$clusters)),Count=table(fit.single$clusters))
            ,file='iCluster.count.txt',sep = '\t',quote = F,row.names = F)

################################

fit.single$clusters
s.inds=match(colnames(gene.exp),paste0(clinical$sample_name,sType))
tm=clinical$os[s.inds]
ev=clinical$os_status[s.inds]
plotKMCox(data.frame(tm,ev,paste0('C',fit.single$clusters)))

tm=clinical$pfs[s.inds]
ev=clinical$pfs_status[s.inds]
plotKMCox(data.frame(tm,ev,paste0('C',fit.single$clusters)))

###############mut################
mut=read.csv(file = mut.file,sep = '\t',check.names = F,stringsAsFactors = F,skip = 5)
head(mut)
mut.genes=unique(mut$SYMBOL)
mtrix.mut=rbind()
sams=c()
for(s in unique(mut$Tumor_Sample_Barcode)){
  gsm=mut$SYMBOL[which(mut$Tumor_Sample_Barcode==s)]
  mt=rep(0,length(mut.genes))
  mt[match(gsm,mut.genes)]=1
  sams=c(sams,substr(s,1,15))
  mtrix.mut=rbind(mtrix.mut,mt)
}
colnames(mtrix.mut)=mut.genes
row.names(mtrix.mut)=sams
head(mtrix.mut)
mtrix.mut[1:3,1:3]
mut.rate=apply(mtrix.mut,2,mean)
mtrix.mut2=mtrix.mut[,which(mut.rate>0.02)]
dim(mtrix.mut2)

mut.inds=match(colnames(gene.exp),row.names(mtrix.mut))
mut.group=fit.single$clusters[!is.na(mut.inds)]
mtrix.mut=mtrix.mut[mut.inds[!is.na(mut.inds)],]
mut.rate=apply(mtrix.mut,2,mean)
table(mut.group)
topgenes=c()
for(i in 1:6){
  mmt=mtrix.mut[which(mut.group==i),]
  mut.rate=apply(mmt,2,mean)
  tp=mut.rate[order(-mut.rate)][1:10]
  topgenes=c(topgenes,names(tp))
}
topgenes=unique(topgenes)
length(topgenes)
all.rate=cbind()
for(i in 1:6){
  mmt=mtrix.mut[which(mut.group==i),]
  mut.rate=apply(mmt,2,mean)
  all.rate=cbind(all.rate,mut.rate[match(topgenes,names(mut.rate))])
}
colnames(all.rate)=paste0('C',1:6)
library(pheatmap)
pheatmap(all.rate,cluster_cols = T,scale = 'row')
dim(all.rate)


#################################################################
rnaseq.count=read.csv(rnaseq.count.file,sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
rnaseq.count=subset(rnaseq.count,select=colnames(gene.exp))
save(rnaseq.count,file = 'rnaseq.count.RData')
dim(rnaseq.count)
gene.type=read.csv(gene.type.file,sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
head(gene.type)
gene.type=gene.type[match(gsub('\\..*','',row.names(rnaseq.count)),row.names(gene.type)),]
pcg.rnaseq.count=rnaseq.count[which(gene.type$TYPE=='protein_coding'),]
lnc.rnaseq.count=rnaseq.count[which(gene.type$TYPE=='lncRNA'),]
cl=sort(unique(fit.single$clusters))


library(DESeq2)
lnc.all=rbind()
pcg.all=rbind()

for(g in cl){
  group=ifelse(fit.single$clusters==g,'Clus','Non')
  dat=lnc.rnaseq.count[,order(group)]
  group=group[order(group)]
  colData <- data.frame(row.names=colnames(dat), group)
  dds <- DESeqDataSetFromMatrix(dat, colData, design= ~ group)
  dds <- dds[rowSums(counts(dds))/ncol(dds) > 1, ]
  dds <- DESeq(dds)
  res= results(dds)
  save(res,file = paste0('C',g,'_lnc_deseq2.RData'))
  al=cbind(row.names(res),res$log2FoldChange,res$padj,rep(g,nrow(res)))
  lnc.all=rbind(lnc.all,al)
  
  group=ifelse(fit.single$clusters==g,'Clus','Non')
  dat=pcg.rnaseq.count[,order(group)]
  group=group[order(group)]
  colData <- data.frame(row.names=colnames(dat), group)
  dds <- DESeqDataSetFromMatrix(dat, colData, design= ~ group)
  dds <- dds[rowSums(counts(dds))/ncol(dds) > 1, ]
  dds <- DESeq(dds)
  res= results(dds)
  save(res,file = paste0('C',g,'_pcg_deseq2.RData'))
  al=cbind(row.names(res),res$log2FoldChange,res$padj,rep(g,nrow(res)))
  pcg.all=rbind(pcg.all,al)
}
save(lnc.all,file = 'lnc.all.deseq.RData')
save(pcg.all,file = 'pcg.all.deseq.RData')
load('lnc.all.deseq.RData')
load('pcg.all.deseq.RData')

pcg.all.dif=which(abs(as.numeric(pcg.all[,2]))>1&as.numeric(pcg.all[,3])<0.05)
lnc.all.dif=which(abs(as.numeric(lnc.all[,2]))>1&as.numeric(lnc.all[,3])<0.05)

save(pcg.all.dif,file = 'pcg.all.dif.deseq.RData')
save(lnc.all.dif,file = 'lnc.all.dif.deseq.RData')
getwd()

load('pcg.all.dif.deseq.RData')
load('lnc.all.dif.deseq.RData')

table((pcg.all[,4]))

dif.count=cbind()
for(g in 0:5){
  dif=pcg.all[pcg.all.dif,][which(pcg.all[pcg.all.dif,4]==g),]
  pcg.c=c(sum(as.numeric(dif[,2])>0),sum(as.numeric(dif[,2])<0),nrow(dif))
  dif=lnc.all[lnc.all.dif,][which(lnc.all[lnc.all.dif,4]==g),]
  lnc.c=c(sum(as.numeric(dif[,2])>0),sum(as.numeric(dif[,2])<0),nrow(dif))
  dif.count=cbind(dif.count,c(pcg.c,lnc.c))
}
colnames(dif.count)=paste0('C',1:6)
row.names(dif.count)=c('PCG_Down','PCG_Up','PCG_All','Lnc_Down','Lnc_Up','Lnc_All')
write.table(cbind(Type=row.names(dif.count),dif.count),'Dif.count.txt',sep='\t',quote=F,row.names=F)

length(unique(gsub('\\..*','',pcg.all[pcg.all.dif,1])))
length(unique(gsub('\\..*','',lnc.all[lnc.all.dif,1])))

group=ifelse(fit.single$clusters==1,'Clus','Non')

pdf('Figure2A.pdf',width = 8,height = 6)
par(mfrow=c(2,3))
for(g in 1:6){
  dif=lnc.all[lnc.all.dif,][which(lnc.all[lnc.all.dif,4]==g),]
  t1=lnc.all[which(lnc.all[,4]==g),]
  lfc=-as.numeric(t1[,2])
  fdr=-log10(as.numeric(t1[,3]))
  col=ifelse(fdr> -log10(0.05),ifelse(lfc>1,'red',ifelse(lfc< -1,'green','blue')),'blue')
  plot(lfc,fdr,pch=16,col=col,ylim=c(0,6),xlim=c(-4,4)
       ,xlab='log2(foldchange)',ylab='-log10(FDR)',main=paste0('SubType ',g))    
}
dev.off()

#fcs=apply(lnc.rnaseq.count,1,function(x){
#  fc=mean(x[which(fit.single$clusters==1)])/mean(x[which(fit.single$clusters!=1)])
#  return(fc)
#})

#tmp=lnc.all[which(lnc.all[,4]=='1'),]
#tmp=tmp[which(abs(as.numeric(tmp[,2]))>1&as.numeric(tmp[,3])<0.05),]
#tmp

#tmp1=cbind(log2(fcs),tmp[match(names(fcs),tmp[,1]),2])[which(log2(fcs)>1),]
#tmp1=tmp1[!is.na(tmp1[,2]),]

#library(pheatmap)
#pheatmap(log2(lnc.rnaseq.count[match(row.names(tmp1),row.names(lnc.rnaseq.count))
#                          ,c(which(fit.single$clusters==1),which(fit.single$clusters!=1))]+1)
#         ,scale = 'row'
#         ,cluster_cols = F
#         )
#
#dif=pcg.all[pcg.all.dif,][which(pcg.all[pcg.all.dif,4]==1),]
#c1.diff=pcg.all[which(pcg.all[,4]==1),]
#dif=data.frame(gene=gs[1:300])
#gs=c1.diff[order(abs(as.numeric(c1.diff[,2]))),1]
#geneLis=abs(as.numeric(c1.diff[order(abs(as.numeric(c1.diff[,2]))),2]))
#names(geneLis)=gs
#geneLis = sort(geneLis, decreasing = TRUE)
#library('clusterProfiler')
library('GSEABase')
#gs1=GeneSet(geneIds=as.character(dif[,1]), setName="subset1")
#all.gmt=list()
for(g in 1:6){
  dif=pcg.all[pcg.all.dif,][which(pcg.all[pcg.all.dif,4]==g),]
  gs1=GeneSet(geneIds=as.character(dif[,1]), setName=paste0('SubType ',g))
  gsc <- GeneSetCollection(list(gs1))
  toGmt(gsc, paste0('GSEA/pcg_C',g,'.gmt'))
  t1=pcg.all[which(pcg.all[,4]==g),]
  lfc=abs(as.numeric(t1[,2]))
  write.table(cbind(Gene=t1,lfc=lfc)[order(lfc),],file = paste0('GSEA/pcg_C',g,'.rnk')
              ,sep = '\t',quote = F,row.names = F)
  dif=lnc.all[lnc.all.dif,][which(lnc.all[lnc.all.dif,4]==g),]
  gs1=GeneSet(geneIds=as.character(dif[,1]), setName=paste0('SubType ',g))
  gsc <- GeneSetCollection(list(gs1))
  #all.gmt[g]=gs1
  toGmt(gsc, paste0('GSEA/lnc_C',g,'.gmt'))
  t1=lnc.all[which(lnc.all[,4]==g),]
  lfc=abs(as.numeric(t1[,2]))
  write.table(cbind(Gene=t1,lfc=lfc)[order(lfc),],file = paste0('GSEA/lnc_C',g,'.rnk')
              ,sep = '\t',quote = F,row.names = F)
}

lnc2disease=read.csv(file='lnc2disease.txt',sep='\t',stringsAsFactors = F)
head(lnc2disease)
lnc2disease$ENSGID
comman.disease.dif=intersect(lnc2disease$ENSGID,unique(gsub('\\..*','',lnc.all[lnc.all.dif,1])))
comman.desease.lnc=intersect(lnc2disease$ENSGID,unique(gsub('\\..*','',lnc.all[,1])))
length(unique(gsub('\\..*','',lnc.all[,1])))
library(VennDiagram)
gene.type[match(comman.disease.dif,row.names(gene.type)),1]
lnc2disease

write.table(cbind(ENSG=comman.disease.dif,Symbol=gene.type[match(comman.disease.dif,row.names(gene.type)),1])
            ,file = 'disease.DElncRNA.txt',row.names = F,quote = F,sep = '\t')

vn=venn.diagram(list(diseaselnc=comman.desease.lnc,DElnc=unique(gsub('\\..*','',lnc.all[lnc.all.dif,1])))
             ,filename = NULL,fill=c("red","yellow"),alpha=c(0.5,0.5))
pdf(file = 'venn.pdf',width = 4,height = 4)
grid.draw(vn)
dev.off()

length(unique(gsub('\\..*','',pcg.all[pcg.all.dif,1])))

fisher.test(matrix(c(length(comman.disease.dif)
                     ,length(unique(gsub('\\..*','',lnc.all[lnc.all.dif,1])))-length(comman.disease.dif)
                     ,length(comman.desease.lnc)-length(comman.disease.dif)
                     ,length(unique(gsub('\\..*','',lnc.all[,1])))
                     -length(comman.desease.lnc)-length(unique(gsub('\\..*','',lnc.all[lnc.all.dif,1])))+-length(comman.disease.dif)
                       ),2,2))
install.packages('UpSetR')
library(UpSetR)
all.lncs=unique(lnc.all[,1])
all.ur=cbind()
for(g in 1:6){
  dif=lnc.all[lnc.all.dif,][which(lnc.all[lnc.all.dif,4]==g),1]
  al=rep(0,length(all.lncs))
  al[match(dif,all.lncs)]=1
  all.ur=cbind(all.ur,al)
}
row.names(all.ur)
colnames(all.ur)=paste0('Subset ',1:6)
all.ur=all.ur[apply(all.ur, 1, sum)>0,]
pdf(file = 'upsetR.pdf',width = 8,height = 4)
upset(as.data.frame(all.ur), nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()
##########################WGCNA###################################
all.tpm=read.csv('Figure/Merge_TCGA-OV_TPM.txt',sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
row.names(all.tpm)

all.dif.lnc=unique(lnc.all[lnc.all.dif,1])
all.dif.pcg=unique(pcg.all[pcg.all.dif,1])
all.tpm=all.tpm[match(gsub('\\..*','',c(all.dif.lnc,all.dif.pcg)),row.names(all.tpm)),]
dim(all.tpm)
gtype.group=c(rep('lnc',length(all.dif.lnc)),rep('pcg',length(all.dif.pcg)))
datExpr=as.data.frame(t(all.tpm))
library(WGCNA)
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
abline(h=5000,col='red')
tpm.inds=which(cutree(sampleTree,h=5000)==1)
datExpr=datExpr[tpm.inds,]
dim(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
pow=sft$powerEstimate
net = blockwiseModules(datExpr, power = pow, maxBlockSize = 7000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)
table(net$colors)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors", 
                                    "GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
wgcna.all.lnc.count=rbind()
for(m in unique(mergedColors)){
#  m='blue'
  t.inds=which(mergedColors==m)
  nlnc=sum(gtype.group[t.inds]=='lnc')
  npcg=sum(gtype.group[t.inds]=='pcg')
  p=fisher.test(matrix(c(nlnc,npcg,length(all.dif.lnc)-nlnc,length(all.dif.pcg)-npcg),2,2)
                ,alternative = 'greater')$p.value
  fc=nlnc*length(all.dif.pcg)/(npcg*length(all.dif.lnc))
  wgcna.all.lnc.count=rbind(wgcna.all.lnc.count,c(m,length(t.inds),nlnc,npcg,p,fc))
}
colnames(wgcna.all.lnc.count)=c('Module','All','Lnc','PCG','p.value','fc')

enrich.lnc.module=wgcna.all.lnc.count[which(as.numeric(wgcna.all.lnc.count[,5])<0.05),]

write.table(wgcna.all.lnc.count,file='wgcna.all.lnc.count.txt',sep = '\t',quote = F,row.names = F)

library('hgu133plus2.db')
all.exp.genes=select(hgu133plus2.db,row.names(all.exp)
                     ,keytype = 'PROBEID',columns=c('ENTREZID','SYMBOL','GENENAME','ENSEMBL'))
sig.module=enrich.lnc.module[-1,1]
library(clusterProfiler)
node.link=rbind()
for(s in sig.module){
  #s='green'
  #s=sig.module[4]
  gs=colnames(datExpr)[which(mergedColors==s)]
  genes=select(hgu133plus2.db,gs
                       ,keytype = 'ENSEMBL',columns=c('ENTREZID','SYMBOL','GENENAME','ENSEMBL'))
  #kegg=enrichKEGG(genes[,2],organism = "hsa",pvalueCutoff = 0.5)
  #kegg=setReadable(kegg, 'hgu133plus2.db',keytype = "ENTREZID")
  go.bp=enrichGO(genes[,2],'hgu133plus2.db',ont = 'ALL',pvalueCutoff = 0.5)
  go.bp=setReadable(go.bp, 'hgu133plus2.db',keytype = "ENTREZID")
  #go.cc=enrichGO(genes[,2],'hgu133plus2.db',ont = 'CC',pvalueCutoff = 0.5)
  #go.cc=setReadable(go.cc, 'hgu133plus2.db',keytype = "ENTREZID")
  #go.mf=enrichGO(genes[,2],'hgu133plus2.db',ont = 'MF',pvalueCutoff = 0.5)
  #go.mf=setReadable(go.mf, 'hgu133plus2.db',keytype = "ENTREZID")
  unique(genes[,2])
  go.bp@result=go.bp@result[go.bp@result$pvalue<0.05,]
  
  summary(go.bp)
  
  #pdf(paste0('Modlue_Enrich/',s,'_go.pdf'),width=8,height=6)
  #dotplot(go.bp,showCategory=20,colorBy='pvalue')
  #barplot(go.bp,showCategory=10)
  #dev.off()
  write.table(summary(go.bp),file = paste0('Modlue_Enrich/',s,'_go.txt'),row.names = F,quote = F,sep = '\t')
  #summary(kegg)$Description
  node.link=rbind(node.link,cbind(rep(s,nrow(go.bp)),summary(go.bp)$Description,summary(go.bp)$ONTOLOGY,summary(go.bp)$ID))
  #summary(go.bp)
  #summary(go.cc)
  #summary(go.mf)
}
colnames(node.link)=c('Module','Pathway','ONT','ID')
write.table(node.link,file = 'Modlue_Enrich/all_module.go.txt',row.names = F,quote = F,sep = '\t')
#rbind(node.link,c(''))
write.table(rbind(node.link,
                  cbind(unique(node.link[,1])
                        ,unique(node.link[,1])
                        ,unique(node.link[,1])
                        ,unique(node.link[,1])))
            ,file = 'Modlue_Enrich/all_module.go.node.txt',row.names = F,quote = F,sep = '\t')



cbind(as.character(summary(go.bp)$ONTOLOGY),summary(go.bp)$ONTOLOGY)

rbind(node.link,
cbind(unique(node.link[,1]),unique(node.link[,1]),unique(node.link[,1]),unique(node.link[,1])))


##############################
cnv.result=read.csv('cnv_results/all_thresholded.by_genes.txt',sep = '\t'
                    ,stringsAsFactors = F,row.names = 1,check.names = F)

cnv.del.result=read.csv('cnv_results/del_genes.conf_90.txt',sep = '\t'
                    ,stringsAsFactors = F,check.names = F)

cnv.amp.result=read.csv('cnv_results/amp_genes.conf_90.txt',sep = '\t'
                        ,stringsAsFactors = F,check.names = F)
gene.type=read.csv(gene.type.file,sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
head(gene.type)

head(cnv.result)
dim(cnv.result)
cnv.result.gids=select(hgu133plus2.db,as.character(cnv.result[,1])
             ,keytype = 'ENTREZID',columns=c('ENTREZID','SYMBOL','GENENAME','ENSEMBL'))

cnv.gene.mp=cbind(cnv.result.gids,gene.type[match(cnv.result.gids[,4],row.names(gene.type)),])
head(cnv.gene.mp)
cnv.gene.mp.lnc=cnv.gene.mp[which(cnv.gene.mp[,6]=='lncRNA'),]
head(cnv.gene.mp.lnc)
length(unique(cnv.gene.mp.lnc[,1]))
lncs=unique(cnv.gene.mp.lnc[,1])
lncs=cnv.gene.mp.lnc[match(lncs,cnv.gene.mp.lnc[,1]),]
head(lncs)

lnc.cnv=cnv.result[match(unique(cnv.gene.mp.lnc[,1]),cnv.result[,1]),c(-1,-2)]
dim(lnc.cnv)
dim(lncs)
lnc.locs=t(apply(lncs,1,function(x){
  crs=unlist(strsplit(x[7],':'))
  sns=unlist(strsplit(crs[2],'-'))
  #print(sns[2])
  end=unlist(strsplit(sns[2],'\\+'))
  #print(end)
  if(length(end)>1) end=end[1]
  return(c(crs[1],sns[1],end))
}))
colnames(lnc.locs)=c('Chromosome','chromStart','chromEnd')
lnc.locs.del=data.frame(Chromosome=as.factor(lnc.locs[,1])
                    ,chromStart=as.numeric(lnc.locs[,2])
                    ,chromEnd=as.numeric(lnc.locs[,3])
                    ,Data=as.numeric((apply(lnc.cnv, 1, function(x){
                      return(sum(x==-2)/length(x))
                    }))))
lnc.locs.ins=data.frame(Chromosome=as.factor(lnc.locs[,1])
                        ,chromStart=as.numeric(lnc.locs[,2])
                        ,chromEnd=as.numeric(lnc.locs[,3])
                        ,Data=as.numeric((apply(lnc.cnv, 1, function(x){
                          return(sum(x==2)/length(x))
                        }))))

library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
data(RCircos.Histogram.Data)
head(RCircos.Histogram.Data)
RCircos.Histogram.Data$Data

RCircos.Set.Core.Components(cyto.info=UCSC.HG38.Human.CytoBandIdeogram
                            ,chr.exclude=NULL,tracks.inside=5,tracks.outside=0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()     
RCircos.Histogram.Plot(lnc.locs.del,data.col = 4,track.num = 2,side = 'in',max.value = 1)
RCircos.Histogram.Plot(lnc.locs.ins,data.col = 4,track.num = 4,max.value = 1)


#dim(cnv.gene.mp.lnc)
#pheatmap(cnv.result[match(unique(cnv.gene.mp.lnc[,1]),cnv.result[,1]),c(-1,-2)]
#         ,cluster_rows = F,cluster_cols = F)

#install.packages('RCircos')

max(lnc.locs.del$Data)
max(lnc.locs.ins$Data)


row.names(lnc.cnv)=gsub('\\|.*','',row.names(lnc.cnv))
lnc.amp=rbind()
for(i in 1:ncol(cnv.amp.result)){
  x=cnv.amp.result[,i]
  x=x[!is.na(x)]
  x=x[4:length(x)]
  cm=intersect(row.names(lnc.cnv),x)
  lnc.amp=rbind(lnc.amp,cbind(rep('Amp',length(cm)),rep(colnames(cnv.amp.result)[i],length(cm)),cm))  
}
for(i in 1:ncol(cnv.del.result)){
  x=cnv.del.result[,i]
  x=x[!is.na(x)]
  x=x[4:length(x)]
  cm=intersect(row.names(lnc.cnv),x)
  lnc.amp=rbind(lnc.amp,cbind(rep('Del',length(cm)),rep(colnames(cnv.del.result)[i],length(cm)),cm))  
}
dim(lnc.cnv)
lnc.amp=cbind(lnc.amp,ENSG=row.names(lncs)[match(lnc.amp[,3],row.names(lnc.cnv))])
head(lnc.amp)

comman.cnv.lncs=intersect(lnc.amp[,4],gsub('\\..*','',all.dif.lnc))
lnc.amp[!is.na(match(lnc.amp[,4],comman.cnv.lncs)),]

write.table(lnc.amp,file = 'lnc.CNV.del.amp.txt',sep = '\t',quote = F,row.names = F)
#cnv.lnc.tpm
#match(lnc.amp[,4],row.names(cnv.lnc.tpm))

table(sort(lnc.amp[which(lnc.amp[,1]=='Del'),2]))
table(sort(lnc.amp[which(lnc.amp[,1]=='Amp'),2]))

dim(all.tpm)

intersrow.names(lncs)
###################################PCC
cnv.lnc.tpm=read.csv('Figure/Merge_TCGA-OV_TPM.txt',sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
cnv.lnc.tpm=cnv.lnc.tpm[match(gsub('\\..*','',row.names(lncs)),gsub('\\..*','',row.names(cnv.lnc.tpm))),]

#dim(cnv.lnc.tpm)
smp.cnv.lnc.tpm=intersect(colnames(cnv.lnc.tpm),colnames(lnc.cnv))
cnv.lnc.tpm=cnv.lnc.tpm[,match(smp.cnv.lnc.tpm,colnames(cnv.lnc.tpm))]
cnv.lnc.cnv=lnc.cnv[,match(smp.cnv.lnc.tpm,colnames(lnc.cnv))]
dim(cnv.lnc.tpm)
dim(cnv.lnc.cnv)
all.pcc=c()
for(i in 1:nrow(cnv.lnc.cnv)){
  all.pcc=c(all.pcc,cor(as.numeric(cnv.lnc.cnv[i,]),as.numeric(cnv.lnc.tpm[i,]))[1])
}
all.pcc.rand=c()
rand.cnv.inds=sample(1:nrow(cnv.lnc.cnv),nrow(cnv.lnc.cnv))
rand.lnc.inds=sample(1:nrow(cnv.lnc.cnv),nrow(cnv.lnc.cnv))
for(i in 1:nrow(cnv.lnc.cnv)){
  all.pcc.rand=c(all.pcc.rand,cor(as.numeric(cnv.lnc.cnv[rand.cnv.inds,][i,]),as.numeric(cnv.lnc.tpm[rand.lnc.inds,][i,]))[1])
}

hist(as.numeric(all.pcc.rand[!is.na(all.pcc.rand)]),col='grey',xlim=c(-1,1),border = 'grey',80)
par(new=T)
plot(density(all.pcc.rand[!is.na(all.pcc.rand)]),xlim=c(-1,1),col='black')
par(new=T)
hist(as.numeric(all.pcc[!is.na(all.pcc)]),col=rgb(255, 0, 0, 80, maxColorValue=255)
     ,xlim=c(-1,1),border = NA,80)
par(new=T)
plot(density(all.pcc[!is.na(all.pcc)]),xlim=c(-1,1),col='blue')
t.test(all.pcc.rand,all.pcc)
############Top 30################
write.table(cbind(lncRNA=row.names(lnc.cnv),lnc.cnv)[which(apply(lnc.cnv, 1, function(x){
  return(sum(abs(x)==2)/length(x))
})>0.3),],file='Figure/P30.lnc.cnv.txt',sep='\t',quote=F,row.names=F)
p30.inds=which(apply(lnc.cnv, 1, function(x){
  return(sum(abs(x)==2)/length(x))
})>0.3)

head(cnv.lnc.cnv)
par(mfrow=c(4,8))
ct=rbind()
#library(beeswarm)
for(i in p30.inds){
  #  i=3
  j=i
  x=cnv.lnc.cnv[j,]
  x1=as.numeric(cnv.lnc.tpm[j,which(x==2)])
  x2=as.numeric(cnv.lnc.tpm[j,which(x==0)])
  x2=c(x2,rep(NA,length(x1)-length(x2)))
  boxplot(data.frame(Amplification=x1,Diploid=x2)
          ,col=c('red','green')
          ,outline=F,ylab='TPM'
          ,main=row.names(cnv.lnc.cnv)[j])
  print(t.test(x1,x2)$p.value)
  legend("topright"
         ,legend=paste0('p=',round(t.test(x1,x2)$p.value,3))
  )
}
##############################################


p03.inds=which(apply(lnc.cnv, 1, function(x){
  return(sum(abs(x)==2)/length(x))
})>0.05&all.pcc>0.1)

#which(all.pcc>0.1)
dim(lnc.cnv)

length(p03.inds)
lnc04=intersect(row.names(lncs)[p03.inds],row.names(all.tpm))
row.names(lncs)[p03.inds][match(lnc04,row.names(lncs)[p03.inds])]

length(row.names(all.tpm))
row.names(lncs)[p03.inds][match(lnc04,row.names(lncs)[p03.inds])]

sig.lncs=cbind(apply(lnc.cnv, 1, function(x){
  return(sum(abs(x)==2)/length(x))
})[p03.inds],all.pcc[p03.inds])[match(lnc04,row.names(lncs)[p03.inds]),]
colnames(sig.lncs)=c('CNV frequancy','PCC')

write.table(cbind(lncRNA=row.names(sig.lncs),sig.lncs),file = 'Figure/sig.lncRNA.txt',sep = '\t',quote = F,row.names = F)

lst.cox.lnc.exp=all.tpm[match(lnc04,gsub('\\..*','',row.names(all.tpm))),]
dim(lst.cox.lnc.exp)
cmm.rnaseq.cli.sample=intersect(colnames(lst.cox.lnc.exp),paste0(clinical$sample_name,sType))
lst.cox.lnc.exp=lst.cox.lnc.exp[,match(cmm.rnaseq.cli.sample,colnames(lst.cox.lnc.exp))]

tm=clinical$os[match(cmm.rnaseq.cli.sample,paste0(clinical$sample_name,sType))]
ev=clinical$os_status[match(cmm.rnaseq.cli.sample,paste0(clinical$sample_name,sType))]

rtm=clinical$pfs[match(cmm.rnaseq.cli.sample,paste0(clinical$sample_name,sType))]
rev=clinical$pfs_status[match(cmm.rnaseq.cli.sample,paste0(clinical$sample_name,sType))]

dim(lst.cox.lnc.exp)
length(tm)
pheatmap(sig.lnc.exp,scale = 'row')

#lnc.cox=t(apply(lst.cox.lnc.exp, 1, function(x){
#  coxFun(data.frame(tm,ev,x)[which(x>0),])  
#}))
lnc.cox=t(apply(lst.cox.lnc.exp, 1, function(x){
  inds=which(x>0)
  tm1=tm[inds]
  ev1=ev[inds]
  x=x[inds]
  dat=data.frame(tm1,ev1,ifelse(x>median(x),'H','L'))
  #dat=dat[which(gp=='Q1'|gp=='Q4'),]
  return(c(plotKMCox(dat),1,1,1))
  #coxFun(data.frame(tm,ev,x)[which(x>0),])  
}))

#which(lnc.cox<0.05)


colnames(lnc.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
lnc.cox=lnc.cox[order(lnc.cox[,1]),]
sig.lnc.cox=lnc.cox[which(lnc.cox[,1]<0.05),]
p.adjust(lnc.cox[,1])
write.table(cbind(LncRNA=row.names(lnc.cox),lnc.cox[,1])[order(lnc.cox[,1]),],file = 'Figure/lst.lnc.cox.txt',sep = '\t',row.names = F,quote = F)


lst.de.counts=lnc.all[lnc.all.dif,][!is.na(match(gsub('\\..*','',lnc.all[lnc.all.dif,1]),row.names(lnc.cox[which(lnc.cox[,1]<0.05),]))),]
cts=cbind()
for(i in 1:6){
  ct=rep(0,nrow(sig.lnc.cox))
  ct[!is.na(match(row.names(sig.lnc.cox)
        ,gsub('\\..*','',lst.de.counts[which(lst.de.counts[,4]==i),1])))]=1
  cts=cbind(cts,ct)
}
row.names(cts)=row.names(sig.lnc.cox)
colnames(cts)=paste0('C',1:6)
library(pheatmap)
pheatmap(cts,cluster_rows = T,cluster_cols = F,color=colorRampPalette(c("navy", "white", "firebrick3"))(50))

apply(cts, 1, sum)

write.table(lncs[match(row.names(cts),row.names(lncs)),],file = 'Figure/lst.lnc.cox.symbol.txt'
            ,sep = '\t',quote = F,row.names = F)

cbind(row.names(cts),all.pcc[match(row.names(cts),row.names(lncs))]
,apply(lnc.cnv[match(row.names(cts),row.names(lncs)),],1,function(x){return(sum(abs(x)==2))})/ncol(lnc.cnv))

cts=cts[which(564*0.05<apply(lnc.cnv[match(row.names(cts),row.names(lncs)),],1,function(x){return(sum(abs(x)==2))})),]

barplot(100*apply(lnc.cnv[match(row.names(cts),row.names(lncs)),],1,function(x){return(sum(abs(x)==2))})/ncol(lnc.cnv)
        ,col='blue',ylab='CNV frequancy')
abline(h=5)

for(i in 1:nrow(cnv.lnc.cnv)){
  all.pcc=c(all.pcc,cor(as.numeric(cnv.lnc.cnv[i,]),as.numeric(cnv.lnc.tpm[i,]))[1])
}

par(mfrow=c(2,3))
for(i in match(row.names(cts),row.names(cnv.lnc.tpm))){
  plot(as.numeric(cnv.lnc.cnv[i,]),as.numeric(cnv.lnc.tpm[i,]),pch=16,xlab='CNV',ylab='lncRNA TPM')
}

########################################
gpl570=read.csv('GPL570-55999.txt',sep = '\t',stringsAsFactors = F)
head(gpl570)
gpl570$ENTREZ_GENE_ID
row.names(lncs)
gpl570=gpl570[!is.na(match(gpl570$ENTREZ_GENE_ID,lncs[,1])),]
dim(gpl570)
cnv2lnc.gpl570=cbind(lncs[match(gpl570$ENTREZ_GENE_ID,lncs[,1]),],gpl570$ID)
write.table(cnv2lnc.gpl570,file = 'GPL5702ENSG.txt',sep = '\t',quote = F,row.names = F)

par(mfrow=c(2,2))
probs=rbind()
for(g in row.names(cts)){
  pb=as.character(cnv2lnc.gpl570[which(cnv2lnc.gpl570[,4]==g),10])
  if(length(pb)==0){
    pb='None'
  }
  probs=rbind(probs,cbind(rep(g,length(pb)),pb))
  x=as.numeric(lst.cox.lnc.exp[match(g,gsub('\\..*','',row.names(lst.cox.lnc.exp))),])
  t.inds=which(x>0)
  tm1=tm[t.inds]
  ev1=ev[t.inds]
  x=x[t.inds]
  #gp=ifelse(x<quantile(x)['25%'],'Q1',ifelse(x<quantile(x)['50%'],'Q2',ifelse(x<quantile(x)['75%'],'Q3','Q4')))
  gp=ifelse(x>median(x),'H','L')
  dat=data.frame(tm1,ev1,gp)
  #dat=dat[which(gp=='Q1'|gp=='Q4'),]
  plotKMCox(dat)
}

probs
###########################################
plotCox=function(dat,pt){
  print(dim(dat))
  colnames(dat)=c('time','status','groups')
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  #print(sdf)
  if(pt==1){
  colKm=rainbow(length(sf$strata))
  plot(sf, mark.time = TRUE,col=colKm,xlab=paste("Survival time in day","\np=",round(p,5)),ylab = "Survival probabilities",main="Method Kaplan Meier")
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(',sf$n,')'), col = colKm,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE)
  }
  return(p)
}

plotKMStep=function(x,time,status,n,pt){
  riskScore=x
  fit <- survivalROC(Stime = time, status = status, marker = x, predict.time = n, method = "KM")
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.1=fit
  d<-c()
  for(j in c(1:length(riskScore))){
    if(riskScore[j]>=optimalCutoff){
      d<-append(d, "High risk")
    }else{
      d<-append(d, "Low risk")
    }
  }
  print(length(which(d=='High risk')))
  #if(length(table(d))==1|min(table(d))<10){
  #  return(1)
  #}else{
  dat=data.frame(time,status,d)
  p=plotCox(dat,pt)
  return(p)
  #}
}
rocKM=function(x,time,status){
  #par(mfrow=c(1,2))
  #time=tm
  #status=ev
  riskScore=as.numeric(x)
  ns=seq(100,3650,100)
  ps=c()
#  for(n in ns){
#    p=plotKMStep(x,time,status,n,0)
#    ps=c(ps,p)
#  }
  print(ps)
#  n=ns[which(ps==min(ps))][1]
#  print(n)
  plotKMStep(x,time,status,3650*3,1)
}

boxplot(t(sig.lnc.exp[match(row.names(cts),gsub('\\..*','',row.names(sig.lnc.exp))),])
        ,outline=F)

apply(sig.lnc.exp[match(row.names(cts),gsub('\\..*','',row.names(sig.lnc.exp))),]
      ,1, function(x){return(sum(x==0))})

par(mfrow=c(3,3))
for(g in row.names(cts)){
  #g=row.names(cts)[4]
  #g=row.names(cts)[5]
  x=as.numeric(sig.lnc.exp[match(g,gsub('\\..*','',row.names(sig.lnc.exp))),])
  t.inds=which(x>0)
  tm1=tm[t.inds]
  ev1=ev[t.inds]
  x=x[t.inds]
  gp=ifelse(x<quantile(x)['25%'],'Q1',ifelse(x<quantile(x)['50%'],'Q2',ifelse(x<quantile(x)['75%'],'Q3','Q4')))
  

  
  #plotKMStep(x,tm1,ev1,365*3,1)
  #rocKM(x,tm,ev)
  #if(g=='ENSG00000253733'){
  #  rocKM((x),tm,ev,600)
  #}else{
  #rocKM(x,tm1,ev1)
  #}
  #0.25*length(x)
  #gp=rep('O',length(x))
  #gp[order(x)[1:100]]='L'
  #gp[order(-x)[1:100]]='H'
  #gp=ifelse(x<=quantile(x)['25%'],'L',ifelse(x>quantile(x)['75%'],'H','O'))
  #print(table(gp))
  #gp=ifelse(x<=mean(x),'L','H')
  #plotKMCox(data.frame(tm1,ev1,gp)[which(gp!='Q2'&gp!='Q1'),])
}


which(t(apply(sig.lnc.exp, 1, function(x){
  coxFun(data.frame(rtm,rev,x))  
}))[,1]<0.05)

plotKMCox(data.frame(tm,ev,paste0('C',fit.single$clusters)))
p3.lnc=row.names(lncs)[which(as.numeric((apply(lnc.cnv, 1, function(x){
  return(sum(abs(x)==2)/length(x))
})))>0.04)]

lest.lnc=intersect(gsub('\\..*','',p3.lnc),gsub('\\..*','',all.dif.lnc))



ct[which(ct[,2]>20),]



max(ct[,1])


dim(all.tpm)
all.de.lnc.cnv=intersect(gsub('\\..*','',row.names(lncs)),gsub('\\..*','',row.names(all.tpm)))

dim(lncs)
p3.sig.lnc.exp=all.tpm[match(all.de.lnc.cnv,gsub('\\..*','',row.names(all.tpm))),]
head(p3.sig.lnc.exp)
p3.sig.lnc.exp=p3.sig.lnc.exp[,match(cmm.rnaseq.cli.sample,colnames(p3.sig.lnc.exp))]

names(which(t(apply(p3.sig.lnc.exp, 1, function(x){
  coxFun(data.frame(tm,ev,x))  
}))[,1]<0.05))

intersect(gsub('\\..*','',lnc.amp[,4]),names(which(t(apply(p3.sig.lnc.exp, 1, function(x){
  coxFun(data.frame(tm,ev,x))  
}))[,1]<0.05)))

#abline(v=mean(all.pcc[!is.na(all.pcc)]))

