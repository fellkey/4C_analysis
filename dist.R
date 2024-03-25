rm(list = ls())
options(stringsAsFactors = F)

getwd()
library(GenomicRanges)

color=c("#D9352B","#4286B7","#5CAB58","#FF874A","#B9BB5A","#8679BE")

###############Function_elements
UCSC=c("T-cell_enhancer_FOMTOM5.bed",
       "T-cell_promoter_FOMTOM5.bed",
       "Jurkat_dbSUPER_SE_UCSC.bed")

name=c('Enhancers','Promoters','Super enhancers')

setwd("~/4c/public_data/rmdup_hg19/")
for (i in 1:3) {
  Func=read.table(UCSC[i])
  Func$V2=as.numeric(Func$V2)
  Func$V3=as.numeric(Func$V3)
  Func=na.omit(Func)
  assign(UCSCname[i], 
         GRanges(seqnames = Func$V1,ranges = IRanges(start = Func$V2,end = Func$V3)))
}

########################IS_analysis_41#####################
setwd("~/4c/newnew/IS_analysis/")
IS=read.csv("./ISannotation_41.csv")

IS_hg19 <- GRanges(seqnames = IS$chr,ranges = IRanges(start = IS$IS_hg19,end = IS$IS_hg19))

#within a linear distance of 2Mb from the HIV-1 integration site.
enhancers=data.frame()
for (i in c(1:41)) {
  distToenhancers=abs(get(name[1])[get(name[1])@seqnames==IS$chr[i]]@ranges@start+
                        get(name[1])[get(name[1])@seqnames==IS$chr[i]]@ranges@width/2-IS$IS_hg19[i])
  distToenhancers
  class(distToenhancers)
  distToenhancers=distToenhancers[distToenhancers<=2e6]
  enhancers1=data.frame(sample=IS$sample[i],status=IS$status[i],distToenhancers=distToenhancers)
  enhancers=rbind(enhancers,enhancers1)
}

promoters=data.frame()
for (i in c(1:41)) {
  distTopromoters=abs(get(name[2])[get(name[2])@seqnames==IS$chr[i]]@ranges@start+
                        get(name[2])[get(name[2])@seqnames==IS$chr[i]]@ranges@width/2-IS$IS_hg19[i])
  distTopromoters
  class(distTopromoters)
  distTopromoters=distTopromoters[distTopromoters<=2e6]
  promoters1=data.frame(sample=IS$sample[i],status=IS$status[i],distTopromoters=distTopromoters)
  promoters=rbind(promoters,promoters1)
}

SE=data.frame()
for (i in c(1:41)) {
  distToSE=abs(get(name[3])[get(name[3])@seqnames==IS$chr[i]]@ranges@start+
                 get(name[3])[get(name[3])@seqnames==IS$chr[i]]@ranges@width/2-IS$IS_hg19[i])
  distToSE
  class(distToSE)
  distToSE=distToSE[distToSE<=2e6]
  SE1=data.frame(sample=IS$sample[i],status=IS$status[i],distToSE=distToSE)
  
  SE=rbind(SE,SE1)
}

#################4C_analysis##############
#within a linear distance of 2Mb from the 4C peaks.
df=read.table("sample_hg19_17.txt")
sample=df$sample
chr=df$chr
IS=df$IS

enhancers=data.frame()
for (i in c(1:17)) {
  Func=read.table(UCSC[1])
  Func$V2=as.numeric(Func$V2)
  Func$V3=as.numeric(Func$V3)
  Func=na.omit(Func)
  Func=Func[Func$V1==chr[i] & Func$V2>=IS[i]-2e6 & Func$V3<=IS[i]+2e6,]
  rangess=GRanges(seqnames = Func$V1,ranges = IRanges(start = Func$V2,end = Func$V3))
  peak4c=peakc[peakc$sample==sample[i] & peakc$start>=IS[i]-2e6 & peakc$end<=IS[i]+2e6,]
  distToenhancers=c(NA)
  for (a in 1:nrow(peak4c)) {
    distToenhancers1=min(abs(rangess@ranges@start+
                               rangess@ranges@width/2-
                               (peak4c$start[a]+peak4c$end[a])/2))
    distToenhancers=c(distToenhancers,distToenhancers1)
  }
  distToenhancers
  class(distToenhancers)
  enhancers1=data.frame(sample=sample[i],chr=chr[i],status=status[i],distToenhancers=distToenhancers)
  enhancers=rbind(enhancers,enhancers1)
}
enhancers=na.omit(enhancers)
enhancers=enhancers[enhancers$distToenhancers!="Inf",]
table(enhancers$sample)

promoters=data.frame()
for (i in c(1:17)) {
  Func=read.table(UCSC[2])
  Func$V2=as.numeric(Func$V2)
  Func$V3=as.numeric(Func$V3)
  Func=na.omit(Func)
  Func=Func[Func$V1==chr[i] & Func$V2>=IS[i]-2e6 & Func$V3<=IS[i]+2e6,]
  rangess=GRanges(seqnames = Func$V1,ranges = IRanges(start = Func$V2,end = Func$V3))
  peak4c=peakc[peakc$sample==sample[i] & peakc$start>=IS[i]-2e6 & peakc$end<=IS[i]+2e6,]
  distTopromoters=c(NA)
  for (a in 1:nrow(peak4c)) {
    distTopromoters1=min(abs(rangess@ranges@start+
                               rangess@ranges@width/2-
                               (peak4c$start[a]+peak4c$end[a])/2))
    distTopromoters=c(distTopromoters,distTopromoters1)
  }
  distTopromoters
  class(distTopromoters)
  promoters1=data.frame(sample=sample[i],chr=chr[i],status=status[i],distTopromoters=distTopromoters)
  promoters=rbind(promoters,promoters1)
}
promoters=na.omit(promoters)
promoters=promoters[promoters$distTopromoters!="Inf",]
table(promoters$sample)

SE=data.frame()
table(peakc$sample)
for (i in c(1:17)) {
  Func=read.table(UCSC[3])
  Func$V2=as.numeric(Func$V2)
  Func$V3=as.numeric(Func$V3)
  Func=na.omit(Func)
  Func=Func[Func$V1==chr[i] & Func$V2>=IS[i]-2e6 & Func$V3<=IS[i]+2e6,]
  rangess=GRanges(seqnames = Func$V1,ranges = IRanges(start = Func$V2,end = Func$V3))
  peak4c=peakc[peakc$sample==sample[i] & peakc$chr==chr[i],]
  distToSE=c(NA)
  for (a in 1:nrow(peak4c)) {
    distToSE1=min(abs(rangess@ranges@start+
                        rangess@ranges@width/2-
                        (peak4c$start[a]+peak4c$end[a])/2))
    distToSE=c(distToSE,distToSE1)
  }
  distToSE
  class(distToSE)
  SE1=data.frame(sample=sample[i],chr=chr[i],status=status[i],distToSE=distToSE)
  SE=rbind(SE,SE1)
}
SE=na.omit(SE)
SE=SE[SE$distToSE!="Inf",]