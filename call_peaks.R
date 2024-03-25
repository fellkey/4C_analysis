rm(list = ls())
options(stringsAsFactors = F)

args <- commandArgs(TRUE)
fasta=args[1]
vpchr=args[2]
vppos=as.numeric(args[3])

setwd(paste0("~/4c/",fasta))

library(peakC)
library(GenomicRanges)
pipe4CFunctionsFile <- "~/4c/pipe4C/functions.R"
source(pipe4CFunctionsFile)

mappedReads <- readRDS(paste0(fasta,".",vpchr,".wSize21.unique.rds"))

vpRegion=4e6
zoom <- GRanges( seqnames=vpchr, resize(IRanges(vppos,vppos),width=vpRegion,fix="center") )
reads <- GRanges(seqnames = mappedReads$chr,ranges = IRanges(start = mappedReads$start,end = mappedRe
ads$end),reads=mappedReads$reads,pos=mappedReads$pos)

vpGR <- reads[unique(queryHits(findOverlaps(reads,zoom)))]
peakCDat <- data.frame(pos=vpGR$pos,reads=vpGR$reads)
wSize=21
alphaFDR=0.05
qWd=1.5
qWr=1
minDist=15e3
resPeakC <- suppressWarnings(single.analysis(data=peakCDat,vp.pos=vppos,wSize=wSize,qWd=qWd,qWr=qWr,minDist=minDist))
resPeakC$vpPos <- vppos
resPeakC$vpChr <- vpchr
min.gapwidth=4e3
vpChr <- resPeakC$vpchr
peakRanges <- reduce(IRanges(resPeakC$peak,resPeakC$peak),min.gapwidth=min.gapwidth)
peakGR <- GRanges(seqnames=vpchr,ranges=peakRanges)

write.table(as.data.frame(peakGR)[,1:3], file=paste0(fasta,"_",vpchr,"_peakC_4mb_unique.txt"), row.names = F, col.names = F, quote=F,sep="\t")
