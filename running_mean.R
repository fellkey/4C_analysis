rm(list = ls())
options(stringsAsFactors = F)

library(GenomicRanges)

args <- commandArgs(TRUE)
fasta=args[1]
chr=args[2]

setwd(paste0("~/4c/",fasta))

AlnReads=read.table(paste0(fasta,".",chr,".unique.bedgraph"))
#head(AlnReads)
AlnReads=AlnReads[,c(1,2,3,8,9)]
colnames(AlnReads)=c("chr","start","end","pos","reads")

nReads=1e6
nTop=2
wSize=21
AlnReads$normReads <- 0
sumTop <- sum( -sort( -AlnReads$reads )[ 1:nTop ] ) 
wNorm <- nReads/( sum( AlnReads$reads )-sumTop ) 
AlnReads$normReads <- wNorm*AlnReads$reads
AlnReads$norm4C <- runmean(x=Rle(AlnReads$normReads),k=wSize,endrule="constant")
AlnReads$norm4C
saveRDS(AlnReads,file=paste0(fasta,".",chr,".wSize21.unique.rds"))
write.table(AlnReads, file=paste0(fasta,".",chr,".wSize21.unique.txt"), row.names = F, col.names = F, quote=F,sep="\t")
